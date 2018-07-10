import sys

from qcip_tools import derivatives, derivatives_e

from nachos.core import fancy_output_derivative, shaking


def get_tensor_converter(type_):
    """Get a converter from ``qcip_tools.derivatives.Tensor`` to a more specialized type"""

    def g(obj):
        if obj.representation.order() == 1:
            kw = {
                'dipole': obj.components
            }
        else:
            kw = {
                'input_fields': tuple(
                    derivatives_e.representation_to_field[a] for a in obj.representation.representation()[1:]),
                'frequency': obj.frequency,
                'tensor': obj.components
            }

        t = type_(**kw)

        return t

    return g


class GetPropertyOfTensor:
    """Get a property for a given object (with an eventual conversion)

    :param what: function
    :type what: callback
    :param convert: convert Tensor to a more specialized tensor
    :type convert: callback
    :param explain: explanation of the property
    :type explain: str
    """
    def __init__(self, accessor, converter=lambda x: x, explain=None, **kwargs):
        self.accessor = accessor
        self.converter = converter
        self.explain = explain
        self.kwargs = kwargs

    def execute(self, obj, total=None):
        """Execute callback on obj, after conversion

        :param obj: the obj
        :type obj: qcip_tools.derivatives.Tensor
        :param total: use callback on ``total - obj`` rather than ``obj``
        :return:
        """
        converted = self.converter(obj)

        if total is not None:
            converted.components = total.components - converted.components

        return self.accessor(converted, **self.kwargs)


def component_access(obj, **kwargs):
    """Access to a coordinate of the tensor

    :param obj: tensor
    :type obj: qcip_tools.derivatives.Tensor
    :rtype: float
    """
    component = kwargs.pop('component')

    if obj.representation.order() != len(component):
        raise ValueError('coordinates shape does not match with the tensor')

    return obj.components[component]


def property_access(obj, **kwargs):
    """Access to a property of the tensor

    :param obj: tensor
    :type obj: qcip_tools.derivatives.Tensor
    """

    func = kwargs.pop('function')

    if not hasattr(obj, func):
        raise ValueError('object {} has no function {}'.format(type(obj), func))

    return getattr(obj, func)(**kwargs)


class BadAnalysis(Exception):
    pass


def _fancy_gather(g, order):
    """Print a partial vibrational contribution (like ``F_F_F``) as fancy one.

    :param g: contributions
    :type g: str
    :param order: perturbation order
    :type order: int
    :rtype: str
    """
    x = g.split('_')
    s = ''
    orders = []
    orders_and_numbers = {}
    for d in x:
        o = len(d)
        if o not in orders:
            orders.append(o)
            orders_and_numbers[o] = 1
        else:
            orders_and_numbers[o] += 1

    orders.sort()
    for i in orders:
        s += shaking.ORDER_TO_REPR[i] + \
            ('' if orders_and_numbers[i] == 1 else shaking.FANCY_EXPONENTS[orders_and_numbers[i]])
    return '[{}]{}'.format(s, shaking.FANCY_EXPONENTS[order])


class Analyzer:
    """Analyzer class that outputs property for different tensors contained in data file

    :param datafile: input for derivatives
    :type datafile: qcip_tools.chemistry_files.chemistry_datafile.ChemistryDataFile
    :param vibrational_contributions: eventual vibrational contributions
    :type vibrational_contributions: dict of nachos.core.shaking.VibrationalContributionsData
    """

    def __init__(self, datafile, vibrational_contributions=None):
        self.datafile = datafile
        self.vibrational_contributions = vibrational_contributions if vibrational_contributions else {}

    def analyze(
            self,
            properties,
            out=sys.stdout,
            only=None,
            frequencies_to_show=None,
            inverse_vibs=None,
            group_vibs=False):
        """Display the different values, as requested

        :param properties: different properties
        :type properties: dict
        :param out: output
        :type out: file
        :param only: derivatives for which properties must be shown
        :type only: list of qcip_tools.derivatives.Derivative
        :param frequencies_to_show: only show given frequencies
        :type frequencies_to_show: list
        :param inverse_vibs: for vibrational, show ``value(total - current)`` rather than ``value(current)``
        :type inverse_vibs: bool
        :param group_vibs: for vibrational, group by perturbation order rather than detailed
        :type group_vibs: bool
        """

        if only:
            bases = only
        else:
            bases = list(
                derivatives.Derivative(a, spacial_dof=self.datafile.spacial_dof)
                for a in self.datafile.derivatives.keys())

            bases.extend(
                derivatives.Derivative(a, spacial_dof=self.datafile.spacial_dof)
                for a in self.vibrational_contributions.keys() if a not in bases)

        bases = [b for b in bases if not derivatives.is_geometrical(b)]
        bases.sort(key=lambda x: (x.order(), x.raw_representation().count('D')))

        no_val = ' ' * 6 + '*' * 3 + ' ' * 6 + ' '

        for base in bases:

            b_repr = base.representation()

            if b_repr not in self.datafile.derivatives and b_repr not in self.vibrational_contributions:
                raise BadAnalysis('{} not in datafile'.format(base))

            if base.order() not in properties:
                continue
            to_show = properties[base.order()]

            frequencies = []

            if b_repr in self.datafile.derivatives:
                frequencies.extend(self.datafile.derivatives[b_repr].keys())

            vibrational_contribution_available = None
            how_much = 2
            gathers = []

            if b_repr in self.vibrational_contributions:
                vibrational_contribution_available = self.vibrational_contributions[b_repr]
                frequencies.extend(vibrational_contribution_available.frequencies_pv)
                if group_vibs:
                    gathers = vibrational_contribution_available.sort_per_type_and_order()
                    how_much += sum(sum(len(gathers[u][e]) for e in gathers[u]) for u in gathers) + 3
                else:
                    how_much += len(vibrational_contribution_available.vibrational_contributions) + 4

            frequencies = list(set(frequencies))
            if frequencies_to_show:
                frequencies_to_show = [derivatives_e.convert_frequency_from_string(i) for i in frequencies_to_show]
                frequencies = list(
                    i for i in frequencies if derivatives_e.convert_frequency_from_string(i) in frequencies_to_show)

            if len(frequencies) == 0:
                continue

            frequencies.sort(key=lambda x: derivatives_e.convert_frequency_from_string(x))

            out.write('\n*** {}:\n'.format(fancy_output_derivative(base)))

            for g in to_show:
                out.write('\n{}:\n'.format(g.explain))
                out.write('-' * 16 * how_much + '\n')

                out.write('{:<15} '.format('frequency'))
                out.write('{:<15} '.format('electronic'))

                if vibrational_contribution_available:
                    inv_mark = ' !!' if inverse_vibs else ''

                    if not group_vibs:
                        for c in vibrational_contribution_available.per_type['zpva']:
                            out.write('{:<15} '.format(c.to_string(fancy=True) + inv_mark))

                        out.write('{:<15} '.format('= ZPVA' + inv_mark))
                    else:
                        if 'zpva' in gathers:
                            out.write('{:<15} '.format('= ZPVA' + inv_mark))

                    if not group_vibs:
                        for c in vibrational_contribution_available.per_type['pv']:
                            out.write('{:<15} '.format(c.to_string(fancy=True) + inv_mark))
                    else:
                        for k in gathers['pv'].keys():
                            for i in range(3):
                                if i not in gathers['pv'][k]:
                                    continue
                                out.write('{:<15} '.format(_fancy_gather(k, i) + inv_mark))

                    out.write('{:<15} '.format('= pv' + inv_mark))
                    out.write('{:<15} '.format('= vib'))
                    out.write('{:<15} '.format('= TOTAL'))

                out.write('\n')
                out.write('-' * 16 * how_much + '\n')

                for frequency in frequencies:
                    out.write('{:<15} '.format(frequency))

                    sum_tensor = derivatives.Tensor(base, frequency=frequency)

                    if b_repr in self.datafile.derivatives and frequency in self.datafile.derivatives[b_repr]:
                        out.write('{: .8e} '.format(g.execute(self.datafile.derivatives[b_repr][frequency])))
                        sum_tensor.components += self.datafile.derivatives[b_repr][frequency].components
                    else:
                        out.write(no_val)

                    if vibrational_contribution_available:
                        v = vibrational_contribution_available.vibrational_contributions

                        # compute sum tensor
                        if frequency in vibrational_contribution_available.frequencies_zpva:
                            sum_tensor.components += \
                                vibrational_contribution_available.total_zpva[frequency].components

                        if frequency in vibrational_contribution_available.frequencies_pv:
                            sum_tensor.components += \
                                vibrational_contribution_available.total_pv[frequency].components

                        # show property values:
                        if frequency in vibrational_contribution_available.frequencies_zpva:
                            if not group_vibs:
                                for c in vibrational_contribution_available.per_type['zpva']:
                                    out.write('{: .8e} '.format(
                                        g.execute(
                                            v[c.to_string()][frequency], total=sum_tensor if inverse_vibs else None)))

                            out.write('{: .8e} '.format(
                                g.execute(
                                    vibrational_contribution_available.total_zpva[frequency],
                                    total=sum_tensor if inverse_vibs else None)))
                        else:
                            for i in range(len(vibrational_contribution_available.per_type['zpva']) + 1):
                                out.write(no_val)

                        if frequency in vibrational_contribution_available.frequencies_pv:
                            if not group_vibs:
                                for c in vibrational_contribution_available.per_type['pv']:
                                    out.write('{: .8e} '.format(
                                        g.execute(
                                            v[c.to_string()][frequency],
                                            total=sum_tensor if inverse_vibs else None)))
                            else:
                                for k in gathers['pv'].keys():
                                    for i in range(3):
                                        if i not in gathers['pv'][k]:
                                            continue
                                        sum_tensor_x = derivatives.Tensor(base, frequency=frequency)
                                        for c in gathers['pv'][k][i]:
                                            sum_tensor_x.components += v[c.to_string()][frequency].components
                                        out.write('{: .8e} '.format(
                                            g.execute(
                                                sum_tensor_x,
                                                total=sum_tensor if inverse_vibs else None)))

                            out.write('{: .8e} '.format(
                                g.execute(
                                    vibrational_contribution_available.total_pv[frequency],
                                    total=sum_tensor if inverse_vibs else None)))
                        else:
                            for i in range(len(vibrational_contribution_available.per_type['pv']) + 1):
                                out.write(no_val)

                        if frequency in vibrational_contribution_available.frequencies_zpva:
                            out.write('{: .8e} '.format(
                                g.execute(vibrational_contribution_available.total_vibrational[frequency])))
                        else:
                            out.write(no_val)

                        out.write('{: .8e} '.format(g.execute(sum_tensor)))

                    out.write('\n')

                out.write('-' * 16 * how_much + '\n')
