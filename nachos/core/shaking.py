import itertools
import math
import sys
import h5py

from nachos.core import fancy_output_derivative

from qcip_tools import derivatives_g, derivatives_e, derivatives
from qcip_tools.chemistry_files import chemistry_datafile


def _merge_dict_of_tensors(b, a):
    """Merge two dicts of tensors. Please kept that function internal."""
    for k in b:
        if k not in a:
            continue
        b[k].components += a[k].components


class VibrationalContributionsData:
    """Store the different vibrational contributions for a given derivative

    :param derivative: the derivative
    :type derivative: qcip_tools.derivatives.Derivative
    """

    def __init__(self, derivative, spacial_dof):
        if derivatives.is_geometrical(derivative):
            raise ValueError('vibrational contributions only to electrical derivatives')

        self.derivative = derivative
        self.vibrational_contributions = {}
        self.spacial_dof = spacial_dof

        self.per_type = {'zpva': [], 'pv': []}

        self.total_zpva = {}
        self.frequencies_zpva = []

        self.total_pv = {}
        self.frequencies_pv = []

        self.total_vibrational = {}

    def add_contribution(self, vc, values):
        """add a contribution

        :param vc: vibrational contribution
        :type vc: VibrationalContribution
        :param values: values (per frequencies)
        :type values: dict
        """

        t = 'zpva' if vc.zpva else 'pv'
        if vc in self.per_type[t]:
            raise ValueError('{} already defined for {}'.format(vc.to_string(), self.derivative))

        for k in values:
            if vc.zpva:
                if k not in self.frequencies_zpva:
                    self.frequencies_zpva.append(k)
                    self.total_zpva[k] = derivatives.Tensor(
                        self.derivative, spacial_dof=self.spacial_dof, frequency=k)
                    self.total_vibrational[k] = derivatives.Tensor(
                        self.derivative, spacial_dof=self.spacial_dof, frequency=k)
            else:
                if k not in self.frequencies_pv:
                    self.frequencies_pv.append(k)
                    self.total_pv[k] = derivatives.Tensor(
                        self.derivative, spacial_dof=self.spacial_dof, frequency=k)

        self.per_type[t].append(vc)
        self.vibrational_contributions[vc.to_string()] = values

        if vc.zpva:
            _merge_dict_of_tensors(self.total_zpva, values)
        else:
            _merge_dict_of_tensors(self.total_pv, values)

        _merge_dict_of_tensors(self.total_vibrational, values)

    def write_in_group(self, group):
        """Write in an h5py group

        :param group: the group
        :type group: h5py.Group
        """

        zpva_contribs = \
            group.attrs['zpva_contributions'].split(',') if 'zpva_contributions' in group.attrs else []
        pv_contribs = \
            group.attrs['pv_contributions'].split(',') if 'pv_contributions' in group.attrs else []

        for k in self.vibrational_contributions:
            chemistry_datafile.ChemistryDataFile.write_derivative_in_dataset(
                group, k, self.derivative, self.vibrational_contributions[k])

        zpva_contribs.extend([p.to_string() for p in self.per_type['zpva'] if p.to_string() not in zpva_contribs])
        pv_contribs.extend([p.to_string() for p in self.per_type['pv'] if p.to_string() not in pv_contribs])

        if zpva_contribs:
            group.attrs['zpva_contributions'] = ','.join(zpva_contribs)
        if pv_contribs:
            group.attrs['pv_contributions'] = ','.join(pv_contribs)

    def read_from_group(self, group):
        zpva_contribs = []
        pv_contribs = []

        if 'zpva_contributions' in group.attrs:
            zpva_contribs = group.attrs['zpva_contributions'].split(',')
        if 'pv_contributions' in group.attrs:
            pv_contribs = group.attrs['pv_contributions'].split(',')

        contributions = zpva_contribs + pv_contribs

        for c in contributions:
            if c not in group:
                raise ValueError('{} not found in {}'.format(c, self.derivative))

            vc = VibrationalContribution.from_representation(c)
            values = chemistry_datafile.ChemistryDataFile.read_derivative_from_dataset(group[c], self.derivative)
            self.add_contribution(vc, values)


def save_vibrational_contributions(path, contributions):
    """Save the contributions in an h5file

    :param path: path to the h5 file
    :type path: str
    :param contributions: the vibrational contributions
    :type contributions: dict
    """

    with h5py.File(path, 'a') as f:
        if '/derivatives/' not in f:
            f.create_group('derivatives')

        if 'vibrational_contributions' not in f['derivatives']:
            vib_group = f['derivatives'].create_group('vibrational_contributions')
        else:
            vib_group = f['derivatives']['vibrational_contributions']

        derivatives_available = \
            vib_group.attrs['derivatives_available'].split(',') if 'derivatives_available' in vib_group.attrs else []

        for k in contributions:
            if k not in vib_group:
                d_group = vib_group.create_group(k)
                derivatives_available.append(k)
            else:
                d_group = vib_group[k]

            contributions[k].write_in_group(d_group)

        vib_group.attrs['derivatives_available'] = ','.join(derivatives_available)


def load_vibrational_contributions(path, spacial_dof):
    """Save the contributions in an h5file

    :param path: path to the h5 file
    :type path: str
    :param spacial_dof: dof
    :type spacial_dof: int
    :rtype: dict
    """

    v_contributions = {}

    with h5py.File(path, 'r') as f:
        if '/derivatives/vibrational_contributions' not in f:
            return v_contributions

        g = f['derivatives']['vibrational_contributions']
        if 'derivatives_available' not in g.attrs:
            raise BadShaking('no derivative available field!')

        available = g.attrs['derivatives_available'].split(',')
        for derivative in available:
            if derivative not in g:
                raise BadShaking('{} is not available'.format(derivative))

            try:
                d = derivatives.Derivative(derivative)
            except derivatives.RepresentationError:
                raise BadShaking('wrong derivative {}'.format(derivative))

            try:
                v_contributions[derivative] = VibrationalContributionsData(d, spacial_dof)
                v_contributions[derivative].read_from_group(g[derivative])
            except ValueError as e:
                raise BadShaking('error while reading {}: {}'.format(derivative, str(e)))

    return v_contributions


class BadShaking(Exception):
    pass


class DerivativeNotAvailable(BadShaking):
    def __init__(self, representation, frequency='static'):
        super().__init__('Derivative not available: {} @ {}'.format(representation, frequency))


ORDER_TO_REPR = {1: 'µ', 2: 'α', 3: 'β', 4: 'γ'}
FANCY_EXPONENTS = {0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴'}


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

    def to_string(self, fancy=False):
        if not fancy:
            return '{}__{}_{}'.format('_'.join(d.representation() for d in self.derivatives), self.m, self.n)
        else:
            signs = ''
            orders = []
            orders_and_numbers = {}
            for i in self.derivatives:
                o = i.order()
                if o not in orders:
                    orders.append(o)
                    orders_and_numbers[o] = 0
                orders_and_numbers[o] += 1

            orders.sort()
            for i in orders:
                k = orders_and_numbers[i]
                signs += '{}{}'.format(ORDER_TO_REPR[i], FANCY_EXPONENTS[k] if k > 1 else '')

            return '[{}]{}{}'.format(signs, FANCY_EXPONENTS[self.m], FANCY_EXPONENTS[self.n])

    def __repr__(self):
        return self.to_string()

    @staticmethod
    def from_representation(representation, dof=3):
        """Convert the **not fancy** representation back into class

        :param representation: reprentation of the vibrational contribution (``X_Y_Z__m_n``)
        :type representation: str
        :param dof: dof
        :type dof: int
        :rtype: VibrationalContribution
        """
        if '__' not in representation:
            raise ValueError('no separator in {}'.format(representation))

        info = representation.split('__')
        if len(info) != 2:
            raise ValueError('incorrect representation {} (number of term)'.format(representation))

        derivatives_ = info[0].split('_')
        anharmonicities = info[1].split('_')

        if len(anharmonicities) != 2:
            raise ValueError('incorrect representation {} (number of anharmonicities)'.format(representation))

        try:
            m, n = int(anharmonicities[0]), int(anharmonicities[1])
        except ValueError:
            raise ValueError('incorrect representation {} (conversion of anharmonicities)'.format(representation))

        converted = []

        for i in derivatives_:
            try:
                converted.append(derivatives.Derivative(i))
            except derivatives.RepresentationError:
                raise ValueError('wrong derivative {}'.format(i))

        return VibrationalContribution(converted, m, n, dof)

    def __eq__(self, other):
        return self.to_string() == other.to_string()


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
        self.available_geometrical_derivatives = {}
        self.dynamic_frequencies = []
        self.available_electrical_derivatives = []

        self.dof = 3 * len(self.datafile.molecule)

        self.computable_pv = {
            # polarizability
            2: [
                VibrationalContribution(('F', 'F'), 0, 0, self.dof),
                VibrationalContribution(('F', 'F'), 1, 1, self.dof),
                VibrationalContribution(('F', 'F'), 2, 0, self.dof),
                VibrationalContribution(('F', 'F'), 0, 2, self.dof)
            ],
            # first hyperpolarizability
            3: [
                VibrationalContribution(('F', 'FF'), 0, 0, self.dof),
                VibrationalContribution(('F', 'FF'), 1, 1, self.dof),
                VibrationalContribution(('F', 'FF'), 2, 0, self.dof),
                VibrationalContribution(('F', 'FF'), 0, 2, self.dof),
                VibrationalContribution(('F', 'F', 'F'), 1, 0, self.dof),
                VibrationalContribution(('F', 'F', 'F'), 0, 1, self.dof)
            ],
            # second hyperpolarizability
            4: [
                VibrationalContribution(('FF', 'FF'), 0, 0, self.dof),
                VibrationalContribution(('F', 'FFF'), 0, 0, self.dof),
                VibrationalContribution(('F', 'F', 'FF'), 1, 0, self.dof),
                VibrationalContribution(('F', 'F', 'FF'), 0, 1, self.dof),
            ]
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

                if base_electrical_derivative not in self.available_geometrical_derivatives:
                    self.available_geometrical_derivatives[base_electrical_derivative] = []
                    self.available_electrical_derivatives.append(derivatives.Derivative(
                        base_electrical_derivative if base_electrical_derivative != 'energy' else '',
                        spacial_dof=self.dof))

                if num_N not in self.available_geometrical_derivatives[base_electrical_derivative]:
                    self.available_geometrical_derivatives[base_electrical_derivative].append(num_N)
            else:
                if k not in self.available_geometrical_derivatives:
                    self.available_geometrical_derivatives[k] = []
                    self.available_electrical_derivatives.append(derivatives.Derivative(k, spacial_dof=self.dof))

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

            if base_electrical not in self.available_geometrical_derivatives:
                return False

            if num_N not in self.available_geometrical_derivatives[base_electrical]:
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

        .. math::
            \\begin{align}
            \\lambda^{\pm ij\\ldots}_{xy\\ldots}= [(\\omega_x+\\omega_y+\\ldots)& +(\\omega_i+\\omega_j+\\ldots)]^{-1}
            \\times\\\\
             &[(\\omega_x+\\omega_y+\\ldots)-(\\omega_i+\\omega_j+\\ldots)]^{-1}
            \\end{align}

        :param up: upper argument (optical frequencies: :math:`\\omega_{i}`, ...)
        :type up: float|list|tuple
        :param down: down argument (vibrational frequencies: :math:`\\omega_{x}`, ...)
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
        # return 1, itertools.permutations(shuflable)

    def compute_zpva(self, vc, derivative, frequencies):
        """Compute a ZPVA contribution to a given derivative. It does not uses _create_tensors() since it is possible
        to go along without permutations (and therefore, work with the all tensor as one).

        :param vc: what to compute
        :type vc: VibrationalContribution
        :param derivative: representation
        :type derivative: qcip_tools.derivatives.Derivative
        :param frequencies: list of frequencies
        :type frequencies: list
        :rtype: dict
        """

        if derivatives.is_geometrical(derivative):
            raise BadShaking('cannot compute vibrational contribution of a geometrical derivative')

        if not vc.zpva or vc.perturbation_order > 1:
            raise BadShaking('cannot compute {}'.format(vc))

        return getattr(self, '_compute_zpva_{}{}'.format(vc.m, vc.n))(derivative, frequencies)

    def compute_pv(self, vc, derivative, frequencies, limit_anharmonicity_usage=True):
        """Compute a pure vibrational contribution.

        .. note::

            Expect the callback function to be ``'_compute_' + what + '_component'``,
            and kwargs to looks like ``'t_' + repr``.

        :param vc: what to compute
        :type vc: VibrationalContribution
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

        if derivative.order() != vc.order:
            raise BadShaking('{} does not match for {}'.format(derivative.representation(), vc.to_string()))

        if not self.check_availability(vc, limit_anharmonicity_usage):
            raise BadShaking('unable to compute {}, some derivatives are missing!'.format(vc.to_string()))

        kwargs = {}
        for n in vc.derivatives_needed(limit_anharmonicity_usage=limit_anharmonicity_usage):
            r = n.representation()
            kwargs['t_' + r.lower()] = self.get_tensor(r)

        return self._create_tensors(
            derivative, frequencies, '_compute_{}_component'.format(vc.to_string()), **kwargs)

    def shake(self, only=None, frequencies=None, out=sys.stdout, verbosity_level=0,
              limit_anharmonicity_usage=True):
        """Compute the vibrational contributions

        :param only: restrict to the vibrational contribution to certain derivatives
        :type only: list|tuple of qcip_tools.derivatives.Derivative, int
        :param frequencies: frequencies (if not available, ZPVA will not be computed for those ones)
        :type frequencies: list
        :param out: output if information is needed to be outputed
        :type out: file
        :param verbosity_level: how far should we print information
        :type verbosity_level: int
        :param limit_anharmonicity_usage: limit the usage of anharmonicity to first order
        :type limit_anharmonicity_usage: bool
        :rtype: dict
        """

        # select bases:
        if not only:
            bases = [(a, 2) for a in self.available_electrical_derivatives if a.order() > 1]
        else:
            bases = []
            for i, max_level in only:
                if i not in self.available_electrical_derivatives:
                    full_F = 'F' * i.order()
                    if full_F not in self.available_electrical_derivatives:
                        raise BadShaking('it is impossible to compute ZPVA or pv contribution to {}'.format(i))

                bases.append((i, max_level))

        bases.sort(key=lambda x: (x[0].order(), x[0].raw_representation().count('D')))

        # select frequencies:
        frequencies_for_all = []
        frequencies_for_pv_only = []

        if not frequencies:
            frequencies_for_all = frequencies_for_pv_only = self.dynamic_frequencies.copy()
        else:
            for f in frequencies:
                if f in self.dynamic_frequencies:
                    frequencies_for_all.append(f)
                frequencies_for_pv_only.append(f)

        vibrational_contributions = {}

        # compute:
        for base, max_level in bases:
            b_repr = base.representation()
            is_dynamic = 'D' in b_repr
            computed_pv, computed_ZPVA = False, False

            if verbosity_level >= 1:
                out.write('\n**** Computing vibrational contributions of {}:\n'.format(fancy_output_derivative(base)))

            freqs_zpva = frequencies_for_all if is_dynamic else ['static']
            freqs_pv = frequencies_for_pv_only if is_dynamic else ['static']

            c = VibrationalContributionsData(base, self.dof)

            to_compute = [VibrationalContribution((base,), 1, 0), VibrationalContribution((base,), 0, 1)]
            if base.order() in self.computable_pv:
                to_compute += self.computable_pv[base.order()]

            for vc in to_compute:
                if vc.perturbation_order > max_level:
                    Shaker.display_message(
                        '{} disabled by request, skipping'.format(vc.to_string(fancy=True)),
                        out,
                        verbosity_level)
                    continue

                if self.check_availability(vc, limit_anharmonicity_usage):
                    Shaker.display_message(
                        'computing {}{}'.format(
                            vc.to_string(fancy=True),
                            ' for all frequencies' if is_dynamic else ''),
                        out,
                        verbosity_level)

                    if vc.zpva:
                        t = self.compute_zpva(vc, base, freqs_zpva)
                        computed_ZPVA = True
                        Shaker.output_tensors(base, vc, t, freqs_zpva, out, verbosity_level)
                    else:
                        t = self.compute_pv(vc, base, freqs_pv, limit_anharmonicity_usage)
                        computed_pv = True
                        Shaker.output_tensors(base, vc, t, freqs_pv, out, verbosity_level)

                    c.add_contribution(vc, t)
                else:
                    Shaker.display_message('unable to compute {}, skipping'.format(vc.to_string(fancy=True)))

            # total stuffs
            if verbosity_level >= 2:
                if computed_ZPVA:
                    out.write('\n*** Total ZPVA contribution to {}:\n'.format(fancy_output_derivative(base)))
                    Shaker.output_tensors(base, None, c.total_zpva, freqs_zpva, out, verbosity_level, 'ZPVA')
                if computed_pv:
                    out.write('\n*** Total pv contribution to {}:\n'.format(fancy_output_derivative(base)))
                    Shaker.output_tensors(base, None, c.total_pv, freqs_pv, out, verbosity_level, 'pv')

            if verbosity_level >= 1 and computed_ZPVA and computed_pv:
                out.write(
                    '\n*** Total vibrational contribution (ZPVA+pv) to {}:\n'.format(fancy_output_derivative(base)))
                Shaker.output_tensors(
                    base, None, c.total_vibrational, freqs_zpva, out, verbosity_level, 'total vibrational contribution')

            vibrational_contributions[b_repr] = c

        return vibrational_contributions

    @staticmethod
    def display_message(message, out=sys.stdout, verbosity_level=0):
        """Output a message, if requested

        :param message: the message
        :type message: str
        :param out: output if information is needed to be outputed
        :type out: file
        :param verbosity_level: how far should we print information
        :type verbosity_level: int
        """

        if verbosity_level >= 1:
            out.write('{}(! {})\n'.format('\n' if verbosity_level > 2 else '', message))
            out.flush()

    @staticmethod
    def output_tensors(base, vc, tensors, frequencies, out=sys.stdout, verbosity_level=0, what=''):
        """Output a bit of information if requested

        :param base: base electrical derivative for which the vibrational contribution is computed
        :type base: qcip_tools.derivatives.Derivative
        :param vc: vibrational contribution
        :type vc: VibrationalContribution|None
        :param tensors: list of tensors
        :type tensors: dict
        :param frequencies: frequencies
        :type frequencies: list
        :param out: output if information is needed to be outputed
        :type out: file
        :param verbosity_level: how far should we print information
        :type verbosity_level: int
        :param what: used when no vc is given
        :type what: str
        """

        if (verbosity_level >= 1 and vc is None) or verbosity_level >= 3:
            for freq in frequencies:
                if freq not in tensors:
                    raise ValueError('{} not in tensors?'.format(freq))
                if vc is not None:
                    out.write('\n** Computed {} for {}:\n'.format(
                        vc.to_string(fancy=True),
                        fancy_output_derivative(base, freq)))
                else:
                    out.write('\n** {} for {}:\n'.format(what, fancy_output_derivative(base, freq)))

                out.write(tensors[freq].to_string(threshold=1e-5))

    # --------------------------------------------
    # BELOW, COMPUTATION OF ALL THE CONTRIBUTIONS:
    #           (so, internal stuffs)
    # --------------------------------------------

    def _create_tensors(self, derivative, frequencies, callback, **kwargs):
        """Create a list of tensor, by taking advantage of the fact that it may (?) be easier to compute the same
        contribution for many frequencies.

        The ``callback`` function must be a function of this class, and receive:

        + `input_fields`` as first argument,
        + then ``frequencies``` (as a list of float, sorted),
        + and finally ``**kwargs``.

        .. note::

            + It is more efficient to compute the static version separately (because of permutations)
            + Permutation of the full tensor is not used: the first coordinate is separated from the rest.


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

        input_fields = [derivatives_e.representation_to_field[x] for x in derivative.representation()[1:]]

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
            for j in derivative.inverse_smart_iterator(i):
                for frequency, value in v.items():
                    initial_frequency = frequencies_mapping[frequency]
                    tensors[initial_frequency].components[j] = value

        return tensors

    def _compute_zpva_10(self, derivative, frequencies):
        """Compute the ZPVA contribution from electrical anharmonicity:

        .. math::

            [p]^{1,0} = \\frac{1}{4}\\,\\sum_a \\tdiff{^2p}{Q_a^2}\\,\\omega_a^{-1}

        :param derivative: the derivative
        :type derivative: qcip_tools.derivatives.Derivative
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :rtype: qcip_tools.derivatives.Tensor
        """
        tensors = {}
        b_repr = derivative.representation()

        for frequency in frequencies:
            t_nnx = self.get_tensor('NN' + b_repr, frequency)
            t = derivatives.Tensor(representation=derivative, frequency=frequency)

            for a in self.mwh.included_modes:
                t.components += 1 / self.mwh.frequencies[a] * t_nnx[a, a]

            t.components *= 1 / 4
            tensors[frequency] = t

        return tensors

    def _compute_zpva_01(self, derivative, frequencies):
        """Compute the ZPVA contribution from mechanical anharmonicity:

        .. math::

            [p]^{0,1} = -\\frac{1}{4}\\,\\sum_{ab} F_{abb}\\,\\tdiff{p}{Q_a}\\,\\omega_a^{-2}\\,\\omega_b^{-1}

        :param derivative: the derivative
        :type derivative: qcip_tools.derivatives.Derivative
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :rtype: qcip_tools.derivatives.Tensor
        """

        tensors = {}
        b_repr = derivative.representation()

        for frequency in frequencies:
            t_nx = self.get_tensor('N' + b_repr, frequency)
            t_nnn = self.get_tensor('NNN')
            t = derivatives.Tensor(representation=derivative, frequency=frequency)

            for a in self.mwh.included_modes:
                ccc = 0
                for b in self.mwh.included_modes:
                    ccc += t_nnn[b, b, a] / self.mwh.frequencies[b]

                t.components += 1 / (self.mwh.frequencies[a] ** 2) * ccc * t_nx[a]

            t.components *= -1 / 4
            tensors[frequency] = t

        return tensors

    def _compute_F_F__0_0_component(self, coo, input_fields, frequencies, t_nf):
        """Compute a component of the :math:`[\\mu^2]^{0,0}` contribution

        .. math::

            [\\mu^2]^{0,0} = \\frac{1}{2}\\,\\sum_{\\mathcal{P}_{ij}} \\sum_a
            \\tdiff{\\mu_i}{Q_a}\\,\\tdiff{\\mu_j}{Q_a}\\,\\lb{\\sigma}{a}

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nf: ``NF`` components
        :type t_nf: numpy.ndarray
        :rtype: list of float
        """

        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                tmp = t_nf[a, p[0][0]] * t_nf[a, p[1][0]]
                for f in frequencies:
                    values[f] += Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * tmp

        for f in frequencies:
            values[f] *= 1 / 2 * multiplier

        return values

    def _compute_F_F__1_1_component(self, coo, input_fields, frequencies, t_nf, t_nnf, t_nnn):
        """Compute a component of the :math:`[\\mu^2]^{1,1}` contribution

        .. math::

            \\begin{align}
            [\\mu^2]^{1,1} &= -\\frac{1}{4}\\,\\sum_{\\mathcal{P}_{ij}} \\sum_{abc} F_{abc}\\,
            \\tdiff{^2\\mu_i}{Q_a\\partial Q_b}\\,\\tdiff{\\mu_j}{Q_c}\\,
            \\lb{\\sigma}{ab}\\,\\lb{\\sigma}{c}\\,(\\omega_a^{-1}+\\omega_b^{-1})\\\\
            &+ F_{bcc}\\,\\tdiff{^2\\mu_i}{Q_a\\partial Q_b}\\,\\tdiff{\\mu_j}{Q_a}\\,
            \\lb{\\sigma}{a}\\,\\omega_b^{-2}\\,\\omega_c^{-1}
            \\end{align}

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nf: ``NF`` components
        :type t_nf: numpy.ndarray
        :param t_nnf: ``NNF`` components
        :type t_nnf: numpy.ndarray
        :param t_nnn: ``NNN`` components
        :type t_nnn: numpy.ndarray
        :rtype: list of float
        """

        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                for b in self.mwh.included_modes:

                    tmp_ab1 = (1 / self.mwh.frequencies[a] + 1 / self.mwh.frequencies[b]) * t_nnf[a, b, p[0][0]]
                    tmp_ab2 = self.mwh.frequencies[b] ** -2 * t_nnf[a, b, p[0][0]] * t_nf[a, p[1][0]]

                    for c in self.mwh.included_modes:
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

    def _compute_F_F__2_0_component(self, coo, input_fields, frequencies, t_nnf):
        """Compute a component of the :math:`[\\mu^2]^{2,0}` contribution

         .. math::

            \\begin{align}
            [\\mu^2]^{2,0} &= \\frac{1}{8}\\,\\sum_{\\mathcal{P}_{ij}} \\sum_{ab}
            \\tdiff{^2\\mu_i}{Q_a\\partial Q_b}\\tdiff{^2\\mu_j}{Q_a\\partial Q_b}\\,
            \\lb{\\sigma}{ab}\\,(\\omega_a^{-1}+\\omega_b^{-1})
            \\end{align}

         :param coo: coordinates
         :type coo: tuple|list
         :param input_fields: input fields
         :type input_fields: tuple|list
         :param frequencies: the frequencies
         :type frequencies: list of float|str
         :param t_nnf: ``NNF`` components
         :type t_nnf: numpy.ndarray
         :rtype: list of float
         """
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                for b in self.mwh.included_modes:

                    tmp = 1 / self.mwh.frequencies[a] * t_nnf[a, b, p[0][0]] * t_nnf[a, b, p[1][0]]

                    for f in frequencies:
                        values[f] += \
                            Shaker.lambda_(p[0][1] * f, (self.mwh.frequencies[a], self.mwh.frequencies[b])) * tmp

        for f in frequencies:
            values[f] *= 1 / 4 * multiplier

        return values

    def _compute_F_F__0_2_component(self, coo, input_fields, frequencies, t_nf, t_nnn):
        """Compute a component of the :math:`[\\mu^2]^{0,2}` contribution

        .. math::

            \\begin{align}
            [\\mu^2]^{0,2} &= \\frac{1}{8}\\,\\sum_{\\mathcal{P}_{ij}} \\sum_{abcd}
            \\tdiff{\\mu_i}{Q_c}\\tdiff{\\mu_j}{Q_d}\\times\\\\
            &\\left[F_{aab}\\,F_{bcd}\\,\\lb{\\sigma}{c}\\lb{\\sigma}{d}\\,\\omega_b^{-2}
            +2\\,F_{abc}\\,F_{abd}\\,\\lb{\\sigma}{ab}\\,\\lb{\\sigma}{c}\\,\\lb{\\sigma}{d}\\right]
            \\end{align}

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nf: ``NF`` components
        :type t_nf: numpy.ndarray
        :param t_nnn: ``NNN`` components
        :type t_nnn: numpy.ndarray
        :rtype: list of float
        """
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                mult_a = 1 / self.mwh.frequencies[a]
                for b in self.mwh.included_modes:
                    mult_ab = self.mwh.frequencies[b] ** -2

                    for c in self.mwh.included_modes:

                        mult_abc = mult_a * t_nf[c, p[0][0]]
                        mult_abc_1 = mult_abc * t_nnn[a, a, b] * mult_ab
                        mult_abc_2 = mult_abc * t_nnn[a, b, c]

                        for d in self.mwh.included_modes:
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

    def _compute_F_FF__0_0_component(self, coo, input_fields, frequencies, t_nf, t_nff):
        """Compute a component of the :math:`[\\mu\\alpha]^{0,0}` contribution

        .. math::

            \\begin{align}
            [\\mu\\alpha]^{0,0} &= \\frac{1}{2}\\,\\sum_{\\mathcal{P}_{ijk}} \\sum_a
            \\tdiff{\\mu_i}{Q_a}\\,\\tdiff{\\alpha_{jk}}{Q_a}\\,\\lb{\\sigma}{a}
            \\end{align}

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nf: ``NF`` components
        :type t_nf: numpy.ndarray
        :param t_nff: ``NFF`` components
        :type t_nff: numpy.ndarray
        :rtype: list of float
        """
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                tmp = t_nf[a, p[0][0]] * t_nff[a, p[1][0], p[2][0]]
                for f in frequencies:
                    values[f] += Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * tmp

        for f in frequencies:
            values[f] *= 1 / 2 * multiplier

        return values

    def _compute_F_FF__2_0_component(self, coo, input_fields, frequencies, t_nnf, t_nnff):
        """Compute a component of the :math:`[\\mu\\alpha]^{2,0}` contribution

        .. math::

            \\begin{align}
            [\\mu\\alpha]^{2,0} &= \\frac{1}{8}\\,\\sum_{\\mathcal{P}_{ijk}} \\sum_{ab}
            \\tdiff{^2\\mu_i}{Q_a\\partial Q_b}\\tdiff{^2\\alpha_{jk}}{Q_a\\partial Q_b}\\,
            \\lb{\\sigma}{ab}\\,(\\omega_a^{-1}+\\omega_b^{-1})\\\\
            \\end{align}

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nnf: ``NNF`` components
        :type t_nnf: numpy.ndarray
        :param t_nnff: ``NNFF`` components
        :type t_nnff: numpy.ndarray
        :rtype: list of float
        """
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                for b in self.mwh.included_modes:

                    tmp = 1 / self.mwh.frequencies[a] * t_nnf[a, b, p[0][0]] * t_nnff[a, b, p[1][0], p[2][0]]

                    for f in frequencies:
                        values[f] += \
                            Shaker.lambda_(p[0][1] * f, (self.mwh.frequencies[a], self.mwh.frequencies[b])) * tmp

        for f in frequencies:
            values[f] *= 1 / 4 * multiplier

        return values

    def _compute_F_FF__0_2_component(self, coo, input_fields, frequencies, t_nf, t_nff, t_nnn):
        """Compute a component of the :math:`[\\mu\\alpha]^{0,2}` contribution

        .. math::

            \\begin{align}
            [\\mu\\alpha]^{0,2} &= \\frac{1}{8}\\,\\sum_{\\mathcal{P}_{ijk}} \\sum_{abcd}
            \\tdiff{\\mu_i}{Q_c}\\tdiff{\\alpha_{jk}}{Q_d}\\times\\\\
            &\\left[ F_{aab}\\,F_{bcd}\\,\\lb{\\sigma}{c}\\,\\lb{\\sigma}{d}\\,\\omega_b^{-2}
            +2\\,F_{abc}\\,F_{abd}\\,\\lb{\\sigma}{ab}\\,\\lb{\\sigma}{c}\\,\\lb{\\sigma}{d}\\right]
            \\end{align}

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nf: ``NF`` components
        :type t_nf: numpy.ndarray
        :param t_nff: ``NFF`` components
        :type t_nff: numpy.ndarray
        :param t_nnn: ``NNN`` components
        :type t_nnn: numpy.ndarray
        :rtype: list of float
        """
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                mult_a = 1 / self.mwh.frequencies[a]
                for b in self.mwh.included_modes:
                    mult_ab = self.mwh.frequencies[b] ** -2

                    for c in self.mwh.included_modes:

                        mult_abc = mult_a * t_nf[c, p[0][0]]
                        mult_abc_1 = mult_abc * t_nnn[a, a, b] * mult_ab
                        mult_abc_2 = mult_abc * t_nnn[a, b, c]

                        for d in self.mwh.included_modes:
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

    def _compute_F_FF__1_1_component(self, coo, input_fields, frequencies, t_nf, t_nnf, t_nff, t_nnff, t_nnn):
        """Compute a component of the :math:`[\\mu\\alpha]^{1,1}` contribution

        .. math::

            \\begin{align}
            [\\mu\\alpha]^{1,1} &= -\\frac{1}{8}\\,\\sum_{\\mathcal{P}_{ij}} \\sum_{abc}
            F_{abc}\\times\\\\
            &\\left[\\tdiff{^2\\mu_i}{Q_a\\partial Q_b}\\,\\tdiff{\\alpha_{jk}}{Q_c}+
            \\tdiff{^2\\alpha_{jk}}{Q_a\\partial Q_b}\\,\\tdiff{\\mu_i}{Q_c}\\right]\\times\\\\
            &\\lb{\\sigma}{ab}\\,\\lb{\\sigma}{c}\\,(\\omega_a^{-1}+\\omega_b^{-1}) \\\\
            &+ F_{bcc}\\,\\left[\\tdiff{^2\\mu_i}{Q_a\\partial Q_b}\\,\\tdiff{\\alpha_{jk}}{Q_a}
            +\\tdiff{^2\\alpha_{jk}}{Q_a\\partial Q_b}\\,\\tdiff{\\mu_i}{Q_a}\\right]\\times\\\\
            &\\lb{\\sigma}{a}\\,\\omega_b^{-2}\\,\\omega_c^{-1}
            \\end{align}

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nnf: ``NNF`` components
        :type t_nnf: numpy.ndarray
        :param t_nnff: ``NNFF`` components
        :type t_nnff: numpy.ndarray
        :param t_nf: ``NF`` components
        :type t_nf: numpy.ndarray
        :param t_nff: ``NFF`` components
        :type t_nff: numpy.ndarray
        :param t_nnn: ``NNN`` components
        :type t_nnn: numpy.ndarray
        :rtype: list of float
        """
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                for b in self.mwh.included_modes:

                    f_ab1 = (1 / self.mwh.frequencies[a] + 1 / self.mwh.frequencies[b])
                    f_ab2 = self.mwh.frequencies[b] ** -2

                    tmp_ab1 = f_ab1 * t_nnf[a, b, p[0][0]]
                    tmp_ab3 = f_ab1 * t_nnff[a, b, p[1][0], p[2][0]]

                    tmp_ab2 = f_ab2 * \
                        (t_nnf[a, b, p[0][0]] * t_nff[a, p[1][0], p[2][0]] +
                         t_nnff[a, b, p[1][0], p[2][0]] * t_nf[a, p[0][0]])

                    for c in self.mwh.included_modes:
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

    def _compute_F_F_F__1_0_component(self, coo, input_fields, frequencies, t_nf, t_nnf):
        """Compute a component of the :math:`[\\mu^3]^{1,0}` contribution

        .. math::

            \\begin{align}
            [\\mu^3]^{1,0} &= \\frac{1}{2}\\,\\sum_{\\mathcal{P}_{ijk}} \\sum_{ab}
            \\tdiff{\\mu_i}{Q_a}\\,\\tdiff{^2\\mu_j}{Q_a\\partial Q_b}\\,
            \\tdiff{\\mu_k}{Q_b}\\,\\lb{\\sigma}{a}\\,\\lb{2}{b}
            \\end{align}

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nnf: ``NNF`` components
        :type t_nnf: numpy.ndarray
        :param t_nf: ``NF`` components
        :type t_nf: numpy.ndarray
        :rtype: list of float
        """
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                tmp_a = t_nf[a, p[0][0]]

                for b in self.mwh.included_modes:
                    tmp_ab = tmp_a * t_nnf[a, b, p[1][0]] * t_nf[b, p[2][0]]

                    for f in frequencies:
                        values[f] += \
                            Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * \
                            Shaker.lambda_(p[2][1] * f, self.mwh.frequencies[b]) * tmp_ab

        for f in frequencies:
            values[f] *= 1 / 2 * multiplier

        return values

    def _compute_F_F_F__0_1_component(self, coo, input_fields, frequencies, t_nf, t_nnn):
        """Compute a component of the :math:`[\\mu^3]^{0,1}` contribution

        .. math::

            \\begin{align}
            [\\mu^3]^{0,1} &= -\\frac{1}{6}\\,\\sum_{\\mathcal{P}_{ijk}} \\sum_{abc} F_{abc}
            \\tdiff{\\mu_i}{Q_a}\\,\\tdiff{\\mu_j}{Q_b}\\,\\tdiff{\\mu_k}{Q_c}\\,
            \\lb{\\sigma}{a}\\,\\lb{1}{b}\\,\\lb{2}{c}
            \\end{align}

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nf: ``NF`` components
        :type t_nf: numpy.ndarray
        :param t_nnn: ``NNN`` components
        :type t_nnn: numpy.ndarray
        :rtype: list of float
        """
        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                tmp_a = t_nf[a, p[0][0]]

                for b in self.mwh.included_modes:
                    tmp_ab = tmp_a * t_nf[b, p[1][0]]

                    for c in self.mwh.included_modes:

                        tmp_abc = tmp_ab * t_nf[c, p[2][0]] * t_nnn[a, b, c]

                        for f in frequencies:
                            values[f] += \
                                Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * \
                                Shaker.lambda_(p[1][1] * f, self.mwh.frequencies[b]) * \
                                Shaker.lambda_(p[2][1] * f, self.mwh.frequencies[c]) * tmp_abc

        for f in frequencies:
            values[f] *= -1 / 6 * multiplier

        return values

    def _compute_FF_FF__0_0_component(self, coo, input_fields, frequencies, t_nff):
        """Compute a component of the :math:`[\\alpha^2]^{0,0}` contribution

        .. math::

            [\\alpha^2]^{0,0} = \\frac{1}{8}\\,\\sum_{\\mathcal{P}_{ijkl}} \\sum_a
            \\tdiff{\\alpha_{ij}}{Q_a}\\,\\tdiff{\\alpha_{kl}}{Q_a}\\,\\lb{23}{a}

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nff: ``NFF`` components
        :type t_nff: numpy.ndarray
        :rtype: list of float
        """

        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                tmp = t_nff[a, p[0][0], p[1][0]] * t_nff[a, p[2][0], p[3][0]]
                for f in frequencies:
                    values[f] += Shaker.lambda_((p[2][1] * f, p[3][1] * f), self.mwh.frequencies[a]) * tmp

        for f in frequencies:
            values[f] *= 1 / 8 * multiplier

        return values

    def _compute_F_FFF__0_0_component(self, coo, input_fields, frequencies, t_nf, t_nfff):
        """Compute a component of the :math:`[\\mu\\beta]^{0,0}` contribution

        .. math::

            [\\mu\\beta]^{0,0} = \\frac{1}{6}\\,\\sum_{\\mathcal{P}_{ijkl}} \\sum_a
            \\tdiff{\\mu_i}{Q_a}\\,\\tdiff{\\beta_{jkl}}{Q_a}\\,\\lb{\\sigma}{a}

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nf: ``NF`` components
        :type t_nf: numpy.ndarray
        :param t_nfff: ``NFFF`` components
        :type t_nfff: numpy.ndarray
        :rtype: list of float
        """

        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                tmp = t_nf[a, p[0][0]] * t_nfff[a, p[1][0], p[2][0], p[3][0]]
                for f in frequencies:
                    values[f] += Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * tmp

        for f in frequencies:
            values[f] *= 1 / 6 * multiplier

        return values

    def _compute_F_F_FF__1_0_component(self, coo, input_fields, frequencies, t_nf, t_nff, t_nnf, t_nnff):
        """Compute a component of the :math:`[\\mu^2\\alpha]^{1,0}` contribution

        .. math::

            [\\mu^2\\alpha]^{1,0} = \\frac{1}{4}\\,\\sum_{\\mathcal{P}_{ijkl}} \\ldots

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nff: ``NFF`` components
        :type t_nff: numpy.ndarray
        :param t_nfff: ``NFFF`` components
        :type t_nfff: numpy.ndarray
        :param t_nnf: ``NNF`` components
        :type t_nnf: numpy.ndarray
        :param t_nnff: ``NNFF`` components
        :type t_nnff: numpy.ndarray
        :rtype: list of float
        """

        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                tmp_a = t_nf[a, p[0][0]]
                for b in self.mwh.included_modes:
                    tmp_ab_1 = tmp_a * t_nnff[a, b, p[1][0], p[2][0]] * t_nf[b, p[3][0]]
                    tmp_ab_2 = tmp_a * t_nnf[a, b, p[1][0]] * t_nff[b, p[2][0], p[3][0]]
                    for f in frequencies:
                        values[f] += \
                            Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * \
                            Shaker.lambda_(p[3][1] * f, self.mwh.frequencies[b]) * \
                            tmp_ab_1
                        values[f] += 2 * \
                            Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * \
                            Shaker.lambda_((p[2][1] * f, p[3][1] * f), self.mwh.frequencies[b]) * \
                            tmp_ab_2

        for f in frequencies:
            values[f] *= 1 / 4 * multiplier

        return values

    def _compute_F_F_FF__0_1_component(self, coo, input_fields, frequencies, t_nf, t_nff, t_nnn):
        """Compute a component of the :math:`[\\mu^2\\alpha]^{0,1}` contribution

        .. math::

            [\\mu^2\\alpha]^{0,1} = -\\frac{1}{4}\\,\\sum_{\\mathcal{P}_{ijkl}} \\ldots

        :param coo: coordinates
        :type coo: tuple|list
        :param input_fields: input fields
        :type input_fields: tuple|list
        :param frequencies: the frequencies
        :type frequencies: list of float|str
        :param t_nff: ``NFF`` components
        :type t_nff: numpy.ndarray
        :param t_nfff: ``NFFF`` components
        :type t_nfff: numpy.ndarray
        :param t_nnn: ``NNN`` components
        :type t_nnn: numpy.ndarray
        :rtype: list of float
        """

        values = {}

        for f in frequencies:
            values[f] = .0

        multiplier, unique_elemts = self.get_iterator(coo, input_fields)

        for p in unique_elemts:
            for a in self.mwh.included_modes:
                tmp_a = t_nf[a, p[0][0]]
                for b in self.mwh.included_modes:
                    tmp_ab = tmp_a * t_nf[b, p[1][0]]
                    for c in self.mwh.included_modes:
                        tmp_abc = tmp_ab * t_nff[c, p[2][0], p[3][0]] * t_nnn[a, b, c]
                        for f in frequencies:
                            values[f] += \
                                Shaker.lambda_(p[0][1] * f, self.mwh.frequencies[a]) * \
                                Shaker.lambda_(p[1][1] * f, self.mwh.frequencies[b]) * \
                                Shaker.lambda_((p[2][1] * f, p[3][1] * f), self.mwh.frequencies[c]) * \
                                tmp_abc

        for f in frequencies:
            values[f] *= -1 / 4 * multiplier

        return values
