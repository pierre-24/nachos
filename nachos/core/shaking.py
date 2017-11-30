import itertools
import math

from qcip_tools import derivatives_g, derivatives_e, derivatives


class BadShaking(Exception):
    pass


class DerivativeNotAvailable(BadShaking):
    def __init__(self, representation, frequency='static'):
        super().__init__('Derivative not available: {} @ {}'.format(representation, frequency))


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
        self.computable_pv = {}

        self.__make_availability()
        self.__make_computable()

    def __make_availability(self):
        """make the availability table
        """

        for k in self.datafile.derivatives:
            num_G = k.count('G')

            if num_G > 0:  # we don't care about pure geometrical derivatives here, we want projected ones
                continue

            num_N = k.count('N')
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

    def __make_computable(self):
        """Make the list of computable stuffs"""

        self.computable_pv = {
            'mu2_00': (2, derivatives_e.PolarisabilityTensor, ('NF',)),
            'mu2_11': (2, derivatives_e.PolarisabilityTensor, ('NF', 'NNF', 'NNN')),
            'mu2_20': (2, derivatives_e.PolarisabilityTensor, ('NNF',)),
            'mu2_02': (2, derivatives_e.PolarisabilityTensor, ('NF', 'NNN')),
            'mu_alpha_00': (3, derivatives_e.PolarisabilityTensor, ('NF', 'NFF')),
            'mu_alpha_11': (3, derivatives_e.PolarisabilityTensor, ('NF', 'NNF', 'NFF', 'NNFF', 'NNN')),
            'mu_alpha_20': (3, derivatives_e.PolarisabilityTensor, ('NNF', 'NNFF')),
            'mu_alpha_02': (3, derivatives_e.PolarisabilityTensor, ('NF', 'NFF', 'NNN')),
            'mu3_10': (3, derivatives_e.PolarisabilityTensor, ('NF', 'NNF')),
            'mu3_01': (3, derivatives_e.PolarisabilityTensor, ('NF', 'NNN')),
        }

    def check_availability(self, derivatives, is_zpva, m, n, limit_anharmonicity_usage=True):
        """Check what is available for the computation

        For a given electrical derivative X,

        + To compute ZPVA, you need NX and NNN ([]⁰¹) and NNX ([]¹⁰)
        + To compute the (n,m) pv contribution of X = [YZ] (X<Y<Z), you need
            + the n+1 geometrical derivative of X and Y
            + the m derivative of the hessian (but it is tweaked so that it only requires NNN if any).
          But X is only contains "F"!

        :param derivatives: list of derivatives (mu, alpha, beta ...)
        :type derivatives: list|tuple
        :param is_zpva: it is zpva ? (m+n should be 1, then)
        :type is_zpva: bool
        :param m: order of electrical anharmonicity
        :type m: int
        :param n: order of mechanical anharmonicity
        :type n: int
        :param limit_anharmonicity_usage: limit the usage of anharmonicity to first order
            (so it accepts m>1, n>1 even if the corresponding derivative is not available)
        :type limit_anharmonicity_usage: bool
        :rtype: bool
        """

        if is_zpva and len(derivatives) != 1:
            raise ValueError('ZPVA is on a single derivative')

        for d in derivatives:
            r = d.representation()

            if not is_zpva and r.count('D') != 0:
                raise ValueError('pure vibrational contributions only on static quantities')

            if r not in self.available:
                return False

            for i in range(m + 1):
                if limit_anharmonicity_usage and i + 1 > 1:
                    continue
                if i + 1 not in self.available[r]:
                    return False

        if limit_anharmonicity_usage:
            if n > 0 and 3 not in self.available['energy']:
                return False
        else:
            for i in range(n + 1):
                if i + 2 not in self.available['energy']:
                    return False

        return True

    def get_tensor(self, representation, frequency='static'):
        if representation in self.datafile.derivatives:
            if not derivatives.is_electrical(representation):
                return self.datafile.derivatives[representation].components

            if frequency in self.datafile.derivatives[representation]:
                return self.datafile.derivatives[representation][frequency].components

        raise DerivativeNotAvailable(representation, frequency)

    def _create_tensors(self, tensor_class, derivative, frequencies, callback, **kwargs):
        """Create a list of tensor, by taking advantage of the fact that it may (?) be easier to compute the same
        contribution for many frequencies.

        The ``callback`` function must be a function of this class, and receive:

        + `input_fields`` as first argument,
        + then ``frequencies``` (as a list of float, sorted),
        + and finally ``**kwargs``.

        .. note::

            It is probably more efficient to compute the static version separately.

        :param tensor_class: class of the tensor to create (if possible, a child of BaseElectricalDerivativeTensor)
        :type tensor_class: class
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
            tensors[frequency] = tensor_class(input_fields=input_fields, frequency=frequency)

        frequencies_converted.sort()

        for i in derivative.smart_iterator():
            v = getattr(self, callback)(i, input_fields, frequencies_converted, **kwargs)
            for frequency, value in v.items():
                initial_frequency = frequencies_mapping[frequency]
                for j in derivative.inverse_smart_iterator(i):
                    tensors[initial_frequency].components[j] = value

        return tensors

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
        """Compute a ZPVA contribution to a given derivative

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

        input_fields = [1 if x == 'D' else 0 for x in derivative.representation()[1:]]

        tensor_classes = {
            1: derivatives_e.ElectricDipole,
            2: derivatives_e.PolarisabilityTensor,
            3: derivatives_e.FirstHyperpolarisabilityTensor,
            4: derivatives_e.SecondHyperpolarizabilityTensor
        }

        if derivative.order() not in tensor_classes:
            raise BadShaking('cannot compute ZPVA contribution to {}'.format(derivative.representation()))

        return getattr(self, '_compute_zpva_{}'.format(what))(
            derivative, frequencies, tensor_classes[derivative.order()], input_fields)

    def compute_pv(self, what, derivative, frequencies):
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
        :rtype: dict
        """

        if derivatives.is_geometrical(derivative):
            raise BadShaking('cannot compute vibrational contribution of a geometrical derivative')

        if what not in self.computable_pv:
            raise BadShaking('{} is not available'.format(what))

        order, tensor_class, needed = self.computable_pv[what]

        if derivative.order() != order:
            raise BadShaking('{} does not match order {} for {}'.format(derivative.representation(), order, what))

        kwargs = {}
        for n in needed:
            kwargs['t_' + n.lower()] = self.get_tensor(n)

        return self._create_tensors(
            tensor_class,
            derivative,
            frequencies,
            '_compute_{}_components'.format(what),
            **kwargs)

    # --------------------------------------------
    # BELOW, COMPUTATION OF ALL THE CONTRIBUTIONS:
    # --------------------------------------------

    def _compute_zpva_10(self, derivative, frequencies, tensor_class, input_fields):
        tensors = {}
        b_repr = derivative.representation()

        for frequency in frequencies:
            t_nnx = self.get_tensor('NN' + b_repr, frequency)
            t = tensor_class(input_fields=input_fields, frequency=frequency)

            for a in range(self.mwh.trans_plus_rot_dof, self.mwh.dof):
                t.components += 1 / self.mwh.frequencies[a] * t_nnx[a, a]

            t.components *= 1 / 4
            tensors[frequency] = t

        return tensors

    def _compute_zpva_01(self, derivative, frequencies, tensor_class, input_fields):
        tensors = {}
        b_repr = derivative.representation()

        for frequency in frequencies:
            t_nx = self.get_tensor('N' + b_repr, frequency)
            t_nnn = self.get_tensor('NNN')
            t = tensor_class(input_fields=input_fields, frequency=frequency)

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
