import math
import collections

from qcip_tools import derivatives, numerical_differentiation, derivatives_e, derivatives_g

CONFIG = {
    'gaussian': {
        'types': ['G', 'F'],
        'methods': [
            ('HF', {'G': 2, 'F': 3}),
            ('DFT', {'G': 2, 'F': 3}),
            ('MP2', {'G': 2, 'F': 2}),  # ... but only static polarizability
            ('MP3', {'G': 1, 'F': 1}),
            ('MP4D', {'G': 1, 'F': 1}),
            ('MP4DQ', {'G': 1, 'F': 1}),
            ('MP4SDQ', {'G': 1, 'F': 1}),
            ('MP4', {'G': 1, 'F': 1}),
            ('CCSD', {'G': 1, 'F': 1}),
            ('CCSD(T)', {'G': 0, 'F': 0}),
        ],
        'bases': ['energy', 'G', 'GG', 'F', 'FF', 'FD', 'FFF', 'FDF', 'FDD'],
        'default_for_extra_fields': {
            'memory': '1Gb',
            'procs': 1,
            'convergence': 11,
            'cc_convergence': 11,
            'cphf_convergence': 10,
            'max_cycles': 600,
            'extra_keywords': '',
            'extra_sections': [],
            'gen_basis': '',
            'vshift': 1000,
            'XC': ''
        }
    },
    'dalton': {
        'types': ['G'],
        'methods': [
            ('CCS', {'G': 1, 'F': 4}),
            ('CC2', {'G': 1, 'F': 4}),
            ('CCSD', {'G': 1, 'F': 4}),
            ('CC3', {'G': 1, 'F': 4}),
        ],
        'bases': ['energy', 'G', 'F', 'FF', 'FD', 'FFF', 'FDF', 'FDD', 'FFFF', 'FDFF', 'FDDF', 'FDDd', 'FDDD'],
        'default_for_extra_fields': {
            'max_iteration': 2500,
            'threshold': 1e-6,
            'cc_threshold': 1e-11,
            'dal_name': 'ND',
        }
    }
}


def compute_numerical_derivative_of_tensor(
        recipe, basis, derivative_repr, tensor_func, frequency=None, dry_run=False, **kwargs):
    """

    :param recipe: recipe
    :type recipe: nachos.core.files.Recipe
    :param basis: basis of differentiation (representation for the base tensors)
    :type basis: qcip_tools.derivatives.Derivative
    :param derivative_repr: representation of the further derivatives of the tensor
    :type derivative_repr: qcip_tools.derivatives.Derivative
    :param tensor_func: function to access to the tensor
    :param frequency: frequency if electrical derivative
    :type frequency: float|str
    :param dry_run: do not fill the tensor or perform romberg analysis
    :type dry_run: bool
    :param kwargs: args passed to tensor_func
    :type kwargs: dict
    :return: tensor
    :rtype: qcip_tools.derivatives.Tensor, dict
    """

    b_repr = basis.representation()
    d_repr = derivative_repr.representation()

    if 'D' in d_repr:
        raise Exception('derivation with respect to D is impossible')
    if 'N' in d_repr:
        raise Exception('derivation with respect to N is impossible (but projection is!)')
    if 'F' in d_repr and 'G' in d_repr:
        raise Exception('no mixed F and G derivatives')
    if 'D' in basis.representation() and not frequency:
        raise Exception('dynamic electric field require frequency')

    derivative_order = derivative_repr.order()

    if basis.spacial_dof is None:
        basis.spacial_dof = derivative_repr.spacial_dof

    final_derivative = basis.differentiate(derivative_repr.representation())
    tensor = derivatives.Tensor(final_derivative, spacial_dof=basis.spacial_dof, frequency=frequency)
    romberg_triangles = collections.OrderedDict()

    for derivative_coo in derivative_repr.smart_iterator():
        romberg_triangles[derivative_coo] = collections.OrderedDict()
        inv_derivative_coos = list(derivative_repr.inverse_smart_iterator(derivative_coo))
        each = collections.Counter(derivative_coo)
        deriv_coefs = []

        for ce in each:
            order = each[ce]
            deriv_coefs.append((numerical_differentiation.Coefficients(
                order,
                numerical_differentiation.Coefficients.choose_p_for_centered(order),
                ratio=recipe['ratio'], method='C'), ce))

        k_max = recipe['k_max']

        if recipe['accuracy_level'] < 1:
            if derivative_order == 3 and derivative_coo == (0, 1, 2):
                # lower the accuracy for the ijk tensor component
                c1f = numerical_differentiation.Coefficients(1, 1, ratio=recipe['ratio'], method='F')
                c1 = numerical_differentiation.Coefficients(1, 2, ratio=recipe['ratio'], method='C')

                if recipe['accuracy_level'] == 0:
                    deriv_coefs = [(c1, 0), (c1f, 1), (c1f, 2)]
                else:  # T-REX
                    deriv_coefs = [(c1f, 0), (c1f, 1), (c1f, 2)]

            if derivative_order == 3 and len(each) == 1:  # lower the accuracy for iii tensor components
                k_max -= 1

        inverse = False
        if recipe['type'] == 'F':
            # because E(-F) = E0 - µ.F + 1/2*a.F² + ..
            #         µ(-F) = µ0 - a.F + 1/2*b.F² + ...
            if basis.order() == 0:
                inverse = True if derivative_order % 2 == 0 else False
            else:
                inverse = True if derivative_order % 2 == 1 else False

        for basis_coo in basis.smart_iterator():
            values_per_k = [
                numerical_differentiation.compute_derivative_of_function(
                    deriv_coefs,
                    tensor_func,
                    k,
                    recipe['min_field'],
                    derivative_repr.shape()[0],
                    # kwargs:
                    basis=basis, inverse=inverse, component=basis_coo, frequency=frequency, recipe=recipe, **kwargs)
                for k in range(k_max)
            ]

            if not dry_run:
                trig = numerical_differentiation.RombergTriangle(values_per_k, ratio=recipe['ratio'], r=2)
                romberg_triangles[derivative_coo][basis_coo] = trig
                val = trig.find_best_value()

                for inv_derivative_coo in inv_derivative_coos:
                    for inv_basis_coo in basis.inverse_smart_iterator(basis_coo):
                        if b_repr != '':
                            if 'F' in d_repr:
                                tensor.components[inv_basis_coo][inv_derivative_coo] = val[1]
                            else:
                                tensor.components[inv_derivative_coo][inv_basis_coo] = val[1]
                        else:
                            tensor.components[inv_derivative_coo] = val[1]

    if not dry_run:
        return tensor, romberg_triangles


def fancy_output_derivative(derivative, frequency=None):
    """

    :param derivative: derivative to output
    :type derivative: qcip_tools.derivatives.Derivative
    :param frequency: eventual frequency
    :type frequency: str|float
    :rtype: str
    """

    r = 'energy'

    g_part = derivative.raw_representation(exclude=derivatives.ELECTRICAL_DERIVATIVES)
    e_part = derivative.raw_representation(exclude=derivatives.GEOMETRICAL_DERIVATIVES)

    if e_part != '':
        if e_part in derivatives_e.NAMES:
            r = derivatives_e.NAMES[e_part]
        else:
            r = derivatives.ORDERS[len(e_part)] + ' order electrical derivative (of the energy)'

        if g_part != '':
            r = derivatives.ORDERS[len(g_part)] + ' order geometrical derivative of ' + r

        if frequency is not None and 'D' in e_part:
            r = '{} (w={})'.format(r, frequency)

    elif g_part != '':
        if g_part in derivatives_g.NAMES:
            r = derivatives_g.NAMES[g_part]
        else:
            r = derivatives.ORDERS[len(g_part)] + ' order geometrical derivative (of the energy)'

    return r


def fancy_output_component_of_derivative(derivative, component, geometry=None):
    """Get the representation of the component of a given tensor

    :param derivative: derivative
    :type derivative: qcip_tools.derivatives.Derivative
    :param component: component of the tensor
    :type component: list|tuple
    :param geometry: the geometry
    :type geometry: qcip_tools.molecule.Molecule
    :rtype: str
    """
    r = ''

    if len(component) != derivative.order():
        raise ValueError('length of component does not match tensor representation')

    for index, d in enumerate(derivative.representation()):
        if d in derivatives.ELECTRICAL_DERIVATIVES:
            r += derivatives.COORDINATES[component[index]]
        else:
            if r != '':
                r += ','
            if geometry is None:
                r += str(component[index])
            else:
                num = int(math.floor(index / 3))
                coo = derivatives.COORDINATES[index % 3]
                r += '{}{}{}'.format(geometry[num].symbol, num + 1, coo)

    return r
