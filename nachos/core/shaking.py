import itertools
import math

from qcip_tools import derivatives_g, derivatives_e, derivatives


class BadShaking(Exception):
    pass


class DerivativeNotAvailable(BadShaking):
    def __init__(self, representation, frequency='static'):
        super().__init__('Derivative not available: {} @ {}'.format(representation, frequency))


class VibrationalContribution:
    """Represent a vibrational contribution

    :param derivatives_: list of derivatives_
    :type derivatives_: list|tuple
    :param m: m (electrical anharmonicity)
    :type m: int
    :param n: n (mechanical anharmonicity)
    :type n: int
    """
    def __init__(self, derivatives_, m=0, n=0, dof=3):

        if m < 0 or n < 0:
            raise ValueError('wrong order of anharmonicity (should be larger than 0)')

        self.order = 0
        self.derivatives = []

        for d in derivatives_:
            if derivatives.is_geometrical(d):
                raise ValueError('{} is not allowed in a vibrational contribution'.format(d))
            elif type(d) is str:
                dx = derivatives.Derivative(d, spacial_dof=dof)
                self.derivatives.append(dx)
                self.order += dx.order()
            elif type(d) is derivatives.Derivative:
                self.derivatives.append(d)
                self.order += d.order()
            else:
                raise ValueError('unable to create a derivative from {}'.format(d))

        self.m = m
        self.n = n
        self.perturbation_order = m + n
        self.dof = dof

        self.zpva = len(self.derivatives) == 1

        if not self.zpva:
            for d in self.derivatives:
                if 'D' in d.representation():
                    raise ValueError('D derivative {} in pv contribution'.format(d))

    def derivatives_needed(self, limit_anharmonicity_usage=True, dof=None):
        """Get the list of the derivatives needed to compute this contribution

        :param limit_anharmonicity_usage: limit the usage to first order of mechanical and electrical anharmonicity
        :type limit_anharmonicity_usage: bool
        :rtype: list
        """

        needed = []

        dof = dof if dof is not None else self.dof

        for i in range(self.m + (0 if not self.zpva else 1), self.m + 2):
            if i < 1:
                continue

            if limit_anharmonicity_usage and i > 2:
                break

            already_treated = []

            for d in self.derivatives:
                if d.order() in already_treated:
                    continue
                needed.append(d.differentiate('N' * i, spacial_dof=dof))
                already_treated.append((d.order()))

        if self.n > 0:
            for i in range(self.n - 1, self.n + 1):
                if i < 1:
                    continue

                if limit_anharmonicity_usage and i > 1:
                    break

                needed.append(derivatives.Derivative('NN' + 'N' * i, spacial_dof=dof))

        return needed

    def __repr__(self):
        return '[{}]^({},{})'.format('+'.join(d.representation() for d in self.derivatives), self.m, self.n)


class Shaker:
    """Shaker class to compute vibrational contributions (to electrical derivatives)

    :param datafile: input for derivatives
    :type datafile: qcip_tools.chemistry_files.chemistry_datafile.ChemistryDataFile
    """

    def __init__(self, datafile):
        self.datafile = datafile

        if 'GG' not in self.datafile.derivatives:
            raise BadShaking('no hessian in datafile')

        self.mwh = derivatives_g.MassWeightedHessian(datafile.molecule, self.datafile.derivatives['GG'])
        self.available = {}

        self.dynamic_frequencies = []

        dof = 3 * len(self.datafile.molecule)

        self.computable_pv = {
            # polarizability
            'mu2_00': VibrationalContribution(('F', 'F'), 0, 0, dof),
            'mu2_11': VibrationalContribution(('F', 'F'), 1, 1, dof),
            'mu2_20': VibrationalContribution(('F', 'F'), 2, 0, dof),
            'mu2_02': VibrationalContribution(('F', 'F'), 0, 2, dof),
            # first hyperpolarizability
            'mu_alpha_00': VibrationalContribution(('F', 'FF'), 0, 0, dof),
            'mu_alpha_11': VibrationalContribution(('F', 'FF'), 1, 1, dof),
            'mu_alpha_20': VibrationalContribution(('F', 'FF'), 2, 0, dof),
            'mu_alpha_02': VibrationalContribution(('F', 'FF'), 0, 2, dof),
            'mu3_10': VibrationalContribution(('F', 'F', 'F'), 1, 0, dof),
            'mu3_01': VibrationalContribution(('F', 'F', 'F'), 0, 1, dof),
        }

        self.__make_availability()

    def __make_availability(self):
        """make the availability table, and the frequency list
        """

        meet_freq = False

        def check_and_deduce(nw, old):
            to_rm = []
            for index, val in enumerate(old):
                if val not in nw:
                    to_rm.append(index)

            deleted = 0
            if to_rm:
                for i in to_rm:
                    del old[i - deleted]
                    deleted += 1

        for k in self.datafile.derivatives:
            num_G = k.count('G')

            if num_G > 0:  # we don't care about pure geometrical derivatives here, we want projected ones
                continue

            num_N = k.count('N')
            num_D = k.count('D')

            if num_D > 0:
                if not meet_freq:
                    self.dynamic_frequencies = [f for f in self.datafile.derivatives[k].keys()]
                    meet_freq = True
                else:
                    check_and_deduce(self.datafile.derivatives[k].keys(), self.dynamic_frequencies)

            if num_N > 0:
                base_electrical_derivative = k.replace('N', '')

                if base_electrical_derivative == '':
                    base_electrical_derivative = 'energy'

                if base_electrical_derivative not in self.available:
                    self.available[base_electrical_derivative] = []

                if num_N not in self.available[base_electrical_derivative]:
                    self.available[base_electrical_derivative].append(num_N)
            else:
                if k not in self.available:
                    self.available[k] = []

        # sort dynamic frequencies
        self.dynamic_frequencies.sort(key=lambda x: derivatives_e.convert_frequency_from_string(x))

    def check_availability(self, vc, limit_anharmonicity_usage=True):
        """Check what is available for the computation

        :param vc: vibrational contribution
        :type vc: VibrationalContribution
        :param limit_anharmonicity_usage: limit the usage of anharmonicity to first order
            (so it accepts m>1, n>1 even if the corresponding derivative is not available)
        :type limit_anharmonicity_usage: bool
        :rtype: bool
        """

        needed = vc.derivatives_needed(limit_anharmonicity_usage)

        for n in needed:
            r = n.representation()
            num_N = r.count('N')
            base_electrical = r.replace('N', '')
            if base_electrical == '':
                base_electrical = 'energy'

            if base_electrical not in self.available:
                return False

            if num_N not in self.available[base_electrical]:
                return False

        return True

    def get_tensor(self, representation, frequency='static'):
        if representation in self.datafile.derivatives:
            if not derivatives.is_electrical(representation):
                return self.datafile.derivatives[representation].components

            if frequency in self.datafile.derivatives[representation]:
                return self.datafile.derivatives[representation][frequency].components

        raise DerivativeNotAvailable(representation, frequency)

    @staticmethod
    def lambda_(up, down):
        """Compute the lambda quantity found in the papers of Kirtman and Bishop

        :param up: upper argument (optical frequencies)
        :type up: float|list|tuple
        :param down: down argument (vibrational frequencies)
        :type down: float|list|tuple
        :rtype: float
        """
        sd = down if type(down) is float else sum(down)
        su = up if type(up) is float else sum(up)

        return (sd + su) ** -1 * (sd - su) ** -1

    @staticmethod
    def get_iterator(coordinates, input_fields):
        """Get the iteration over all possible permutations

        :param coordinates: coordinates
        :type coordinates: tuple
        :param input_fields: the input fields
        :type input_fields: tuple
        :rtype: (int, list)
        """

        shuflable = [(coordinates[0], -sum(input_fields))]
        for i in range(1, len(coordinates)):
            shuflable.append((coordinates[i], input_fields[i - 1]))

        unique_elements = set(itertools.permutations(shuflable))
        return math.factorial(len(coordinates)) / len(unique_elements), unique_elements

    def compute_zpva(self, what, derivative, frequencies):
        """Compute a ZPVA contribution to a given derivative. It does not uses _create_tensors() since it is possible
        to go along without permutations (and therefore, work with the all tensor as one).

        :param what: what to compute
        :type what: str
        :param derivative: representation
        :type derivative: qcip_tools.derivatives.Derivative
        :param frequencies: list of frequencies
        :type frequencies: list
        :rtype: dict
        """

        if derivatives.is_geometrical(derivative):
            raise BadShaking('cannot compute vibrational contribution of a geometrical derivative')

        if what not in ['10', '01']:
            raise BadShaking('cannot compute zpva_{}'.format(what))

        t_ = {'10': (1, 0), '01': (0, 1)}

        if not self.check_availability(VibrationalContribution((derivative,), *t_[what])):
            raise BadShaking('unable to compute zpva_{}, some derivatives are missing'.format(what))

        return getattr(self, '_compute_zpva_{}'.format(what))(derivative, frequencies)

    def compute_pv(self, what, derivative, frequencies, limit_anharmonicity_usage=True):
        """Compute a pure vibrational contribution.

        .. note::

            Expect the callback function to be ``'_compute_' + what + '_components'``,
            and kwargs to looks like ``'t_' + repr``.

        :param what: what to compute
        :type what: str
        :param derivative: representation
        :type derivative: qcip_tools.derivatives.Derivative
        :param frequencies: list of frequencies
        :type frequencies: list
        :param limit_anharmonicity_usage: limit the usage of anharmonicity to first order
        :type limit_anharmonicity_usage: bool
        :rtype: dict
        """

        if derivatives.is_geometrical(derivative):
            raise BadShaking('cannot compute vibrational contribution of a geometrical derivative')

        if what not in self.computable_pv:
            raise BadShaking('{} is not available'.format(what))

        vc = self.computable_pv[what]

        if derivative.order() != vc.order:
            raise BadShaking('{} does not match order {} for {}'.format(derivative.representation(), vc.order, what))

        if not self.check_availability(vc, limit_anharmonicity_usage):
            raise BadShaking('unable to compute {}, some derivatives are missing!'.format(what))

        kwargs = {}
        for n in vc.derivatives_needed(limit_anharmonicity_usage=limit_anharmonicity_usage):
            r = n.representation()
            kwargs['t_' + r.lower()] = self.get_tensor(r)

        return self._create_tensors(derivative, frequencies, '_compute_{}_components'.format(what), **kwargs)

    # --------------------------------------------
    # BELOW, COMPUTATION OF ALL THE CONTRIBUTIONS:
    # --------------------------------------------

    def _create_tensors(self, derivative, frequencies, callback, **kwargs):
        """Create a list of tensor, by taking advantage of the fact that it may (?) be easier to compute the same
        contribution for many frequencies.

        The ``callback`` function must be a function of this class, and receive:

        + `input_fields`` as first argument,
        + then ``frequencies``` (as a list of float, sorted),
        + and finally ``**kwargs``.

        .. note::

            It is probably more efficient to compute the static version separately.

        :param derivative: the derivative of the tensor for which the contribution should be computed
        :type derivative: qcip_tools.derivatives.Derivative
        :param frequencies: the frequencies
        :type frequencies: list
        :param callback: callback func
        :type callback: str
        :param kwargs: kwargs
        :type kwargs: dict
        :rtype: dict
        """

        if derivative.representation() not in derivatives_e.DERIVATIVES:
            raise BadShaking('I cannot deal with {}'.format(derivative.representation()))

        input_fields = [1 if x == 'D' else 0 for x in derivative.representation()[1:]]

        tensors = {}
        frequencies_converted = []
        frequencies_mapping = {}

        for frequency in frequencies:
            converted_frequency = derivatives_e.convert_frequency_from_string(frequency)
            frequencies_converted.append(converted_frequency)
            frequencies_mapping[converted_frequency] = frequency
            tensors[frequency] = derivatives.Tensor(representation=derivative, frequency=frequency)

        frequencies_converted.sort()

        for i in derivative.smart_iterator():
            v = getattr(self, callback)(i, input_fields, frequencies_converted, **kwargs)
            for frequency, value in v.items():
                initial_frequency = frequencies_mapping[frequency]
                for j in derivative.inverse_smart_iterator(i):
                    tensors[initial_frequency].components[j] = value

        return tensors

    def _compute_zpva_10(self, derivative, frequencies):
        tensors = {}
        b_repr = derivative.representation()

        for frequency in frequencies:
            t_nnx = self.get_tensor('NN' + b_repr, frequency)
            t = derivatives.Tensor(representation=derivative, frequency=frequency)

            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                t.components += 1 / self.mwh.frequencies[a] * t_nnx[a, a]

            t.components *= 1 / 4
            tensors[frequency] = t

        return tensors

    def _compute_zpva_01(self, derivative, frequencies):
        tensors = {}
        b_repr = derivative.representation()

        for frequency in frequencies:
            t_nx = self.get_tensor('N' + b_repr, frequency)
            t_nnn = self.get_tensor('NNN')
            t = derivatives.Tensor(representation=derivative, frequency=frequency)

            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                ccc = 0
                for b in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                    ccc += t_nnn[b, b, a] / self.mwh.frequencies[b]

                t.components += 1 / (self.mwh.frequencies[a] ** 2) * ccc * t_nx[a]

            t.components *= -1 / 4
            tensors[frequency] = t

        return tensors

    def _compute_mu2_00_components(self, coo, input_fields, frequencies, t_nf):
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                tmp = t_nf[a, p[0][0]] * t_nf[a, p[1][0]]
                for f in frequencies:
                    values[f] += Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * tmp

        for f in frequencies:
            values[f] *= 1 / 2 * multiplier

        return values

    def _compute_mu2_11_components(self, coo, input_fields, frequencies, t_nf, t_nnf, t_nnn):
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                for b in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):

                    tmp_ab1 = (1 / self.mwh.frequencies[a] + 1 / self.mwh.frequencies[b]) * t_nnf[a, b, p[0][0]]
                    tmp_ab2 = self.mwh.frequencies[b] ** -2 * t_nnf[a, b, p[0][0]] * t_nf[a, p[1][0]]

                    for c in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                        tmp1 = tmp_ab1 * t_nf[c, p[1][0]] * t_nnn[a, b, c]
                        tmp2 = tmp_ab2 * t_nnn[b, c, c] * self.mwh.frequencies[c] ** -1

                        for f in frequencies:
                            ws = p[0][1] * f
                            fr_1 = Shaker.lambda_(ws, (self.mwh.frequencies[a], self.mwh.frequencies[b])) * \
                                Shaker.lambda_(ws, self.mwh.frequencies[c])

                            values[f] += (fr_1 * tmp1 + Shaker.lambda_(ws, self.mwh.frequencies[a]) * tmp2)

        for f in frequencies:
            values[f] *= - 1 / 4 * multiplier

        return values

    def _compute_mu2_20_components(self, coo, input_fields, frequencies, t_nnf):
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                for b in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):

                    tmp = 1 / self.mwh.frequencies[a] * t_nnf[a, b, p[0][0]] * t_nnf[a, b, p[1][0]]

                    for f in frequencies:
                        values[f] += \
                            Shaker.lambda_(p[0][1] * f, (self.mwh.frequencies[a], self.mwh.frequencies[b])) * tmp

        for f in frequencies:
            values[f] *= 1 / 4 * multiplier

        return values

    def _compute_mu2_02_components(self, coo, input_fields, frequencies, t_nf, t_nnn):
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                mult_a = 1 / self.mwh.frequencies[a]
                for b in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                    mult_ab = self.mwh.frequencies[b] ** -2

                    for c in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):

                        mult_abc = mult_a * t_nf[c, p[0][0]]
                        mult_abc_1 = mult_abc * t_nnn[a, a, b] * mult_ab
                        mult_abc_2 = mult_abc * t_nnn[a, b, c]

                        for d in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                            mult_abcd = t_nf[d, p[1][0]]

                            tmp1 = mult_abcd * mult_abc_1 * t_nnn[b, c, d]
                            tmp2 = mult_abcd * mult_abc_2 * t_nnn[a, b, d] * 2

                            for f in frequencies:
                                ws = p[0][1] * f
                                values[f] -= \
                                    (tmp1 +
                                     Shaker.lambda_(ws, (self.mwh.frequencies[a], self.mwh.frequencies[b])) * tmp2) * \
                                    self.lambda_(ws, self.mwh.frequencies[c]) * \
                                    self.lambda_(ws, self.mwh.frequencies[d])

        for f in frequencies:
            values[f] *= -1 / 8 * multiplier

        return values

    def _compute_mu_alpha_00_components(self, coo, input_fields, frequencies, t_nf, t_nff):
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                tmp = t_nf[a, p[0][0]] * t_nff[a, p[1][0], p[2][0]]
                for f in frequencies:
                    values[f] += Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * tmp

        for f in frequencies:
            values[f] *= 1 / 2 * multiplier

        return values

    def _compute_mu_alpha_20_components(self, coo, input_fields, frequencies, t_nnf, t_nnff):
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                for b in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):

                    tmp = 1 / self.mwh.frequencies[a] * t_nnf[a, b, p[0][0]] * t_nnff[a, b, p[1][0], p[2][0]]

                    for f in frequencies:
                        values[f] += \
                            Shaker.lambda_(p[0][1] * f, (self.mwh.frequencies[a], self.mwh.frequencies[b])) * tmp

        for f in frequencies:
            values[f] *= 1 / 4 * multiplier

        return values

    def _compute_mu_alpha_02_components(self, coo, input_fields, frequencies, t_nf, t_nff, t_nnn):
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                mult_a = 1 / self.mwh.frequencies[a]
                for b in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                    mult_ab = self.mwh.frequencies[b] ** -2

                    for c in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):

                        mult_abc = mult_a * t_nf[c, p[0][0]]
                        mult_abc_1 = mult_abc * t_nnn[a, a, b] * mult_ab
                        mult_abc_2 = mult_abc * t_nnn[a, b, c]

                        for d in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                            mult_abcd = t_nff[d, p[1][0], p[2][0]]

                            tmp1 = mult_abcd * mult_abc_1 * t_nnn[b, c, d]
                            tmp2 = mult_abcd * mult_abc_2 * t_nnn[a, b, d] * 2

                            for f in frequencies:
                                ws = p[0][1] * f
                                values[f] -= \
                                    (tmp1 +
                                     Shaker.lambda_(ws, (self.mwh.frequencies[a], self.mwh.frequencies[b])) * tmp2) * \
                                    self.lambda_(ws, self.mwh.frequencies[c]) * \
                                    self.lambda_(ws, self.mwh.frequencies[d])

        for f in frequencies:
            values[f] *= -1 / 8 * multiplier

        return values

    def _compute_mu_alpha_11_components(self, coo, input_fields, frequencies, t_nf, t_nnf, t_nff, t_nnff, t_nnn):
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                for b in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):

                    f_ab1 = (1 / self.mwh.frequencies[a] + 1 / self.mwh.frequencies[b])
                    f_ab2 = self.mwh.frequencies[b] ** -2

                    tmp_ab1 = f_ab1 * t_nnf[a, b, p[0][0]]
                    tmp_ab3 = f_ab1 * t_nnff[a, b, p[1][0], p[2][0]]

                    tmp_ab2 = f_ab2 * \
                        (t_nnf[a, b, p[0][0]] * t_nff[a, p[1][0], p[2][0]] +
                         t_nnff[a, b, p[1][0], p[2][0]] * t_nf[a, p[0][0]])

                    for c in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                        tmp1 = (tmp_ab1 * t_nff[c, p[1][0], p[2][0]] + tmp_ab3 * t_nf[c, p[0][0]]) * t_nnn[a, b, c]
                        tmp2 = tmp_ab2 * t_nnn[b, c, c] * self.mwh.frequencies[c] ** -1

                        for f in frequencies:
                            ws = p[0][1] * f
                            fr_1 = Shaker.lambda_(ws, (self.mwh.frequencies[a], self.mwh.frequencies[b])) * \
                                Shaker.lambda_(ws, self.mwh.frequencies[c])

                            values[f] += (fr_1 * tmp1 + Shaker.lambda_(ws, self.mwh.frequencies[a]) * tmp2)

        for f in frequencies:
            values[f] *= - 1 / 8 * multiplier

        return values

    def _compute_mu3_10_components(self, coo, input_fields, frequencies, t_nf, t_nnf):
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                tmp_a = t_nf[a, p[0][0]]

                for b in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                    tmp_ab = tmp_a * t_nnf[a, b, p[1][0]] * t_nf[b, p[2][0]]

                    for f in frequencies:
                        values[f] += \
                            Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * \
                            Shaker.lambda_(p[2][1] * f, self.mwh.frequencies[b]) * tmp_ab

        for f in frequencies:
            values[f] *= 1 / 2 * multiplier

        return values

    def _compute_mu3_01_components(self, coo, input_fields, frequencies, t_nf, t_nnn):
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                tmp_a = t_nf[a, p[0][0]]

                for b in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                    tmp_ab = tmp_a * t_nf[b, p[1][0]]

                    for c in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):

                        tmp_abc = tmp_ab * t_nf[c, p[2][0]] * t_nnn[a, b, c]

                        for f in frequencies:
                            values[f] += \
                                Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * \
                                Shaker.lambda_(p[1][1] * f, self.mwh.frequencies[b]) * \
                                Shaker.lambda_(p[2][1] * f, self.mwh.frequencies[c]) * tmp_abc

        for f in frequencies:
            values[f] *= -1 / 6 * multiplier

        return values
