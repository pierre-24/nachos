import math

from qcip_tools import derivatives, derivatives_e, derivatives_g

CONFIG = {
    'gaussian': {
        'types': ['G', 'F'],
        'methods': [
            ('HF', {'G': 2, 'F': 3}),
            ('DFT', {'G': 2, 'F': 3}),
            ('MP2', {'G': 2, 'F': 2}),  # ... but only static polarizability
            ('SCS-MP2', {'G': 0, 'F': 0}),  # special MP2-SCS, from energy
            ('MP3', {'G': 1, 'F': 1}),
            ('MP4D', {'G': 1, 'F': 1}),
            ('MP4DQ', {'G': 1, 'F': 1}),
            ('MP4SDQ', {'G': 1, 'F': 1}),
            ('MP4', {'G': 1, 'F': 1}),
            ('CCSD', {'G': 1, 'F': 0}),
            ('CCSD(T)', {'G': 0, 'F': 0}),
        ],
        'bases': ['energy', 'G', 'GG', 'F', 'FF', 'dD', 'FFF', 'dDF', 'XDD'],
        'default_for_extra_fields': {
            'memory': '1Gb',
            'procs': 1,
            'convergence': 11,
            'cc_convergence': 11,
            'cphf_convergence': 10,
            'max_cycles': 600,
            'extra_keywords': '',
            'extra_sections': '',
            'extra_sections_before': '',
            'gen_basis': '',
            'vshift': 1000,
            'XC': '',
            'use_full': True,
        }
    },
    'dalton': {
        'types': ['G'],
        'methods': [
            ('HF', {'G': 2, 'F': 4}),
            ('DFT', {'G': 2, 'F': 4}),  # with some XC functionals only :o
            ('CC', {'G': 1, 'F': 4}),
        ],
        'bases': [
            'energy', 'G', 'GG', 'F', 'FF', 'dD', 'FFF', 'dDF', 'XDD', 'FFFF', 'dFFD', 'dDFF', 'XDDF', 'dDDd', 'XDDD'],
        'default_for_extra_fields': {
            'threshold': 1e-11,
            'cc_threshold': 1e-11,
            'response_threshold': 1e-10,
            'response_max_it': 2500,
            'response_max_ito': 10,
            'response_dim_reduced_space': 2500,
            'dal_name': 'ND',
            'CC': '',
            'XC': '',
            'split_level_3': 1,
            'split_level_4': 1,
            'merge_level_3': 0,
            'merge_level_4': 0
        }
    },

    'qchem': {
        'types': ['F', 'G'],
        'methods': [
            # Only works with CCMAN2 (!)
            # See http://www.q-chem.com/qchem-website/manual/qchem51_manual/sect-ccmeth.html
            # for a list of available CCMAN2 methods
            ('MP2', {'G': 0, 'F': 0}),
            ('MP3', {'G': 0, 'F': 0}),
            ('QCISD', {'G': 0, 'F': 0}),
            ('QCISD(T)', {'G': 0, 'F': 0}),
            ('CCD', {'G': 0, 'F': 0}),
            ('CCSD', {'G': 0, 'F': 0}),
            ('CCSD(T)', {'G': 0, 'F': 0})
        ],
        'bases': ['energy'],
        'default_for_extra_fields': {
            'convergence': 11,
            'cc_convergence': 10,
            'max_cycles': 600,
            'memory_static': 2000,
            'memory_cc': 2000
        }
    }
}

GAUSSIAN_DOUBLE_HYBRIDS = [
    'B2PLYP', 'B2PLYPD', 'B2PLYPD3',
    'mPW2PLYP', 'mPW2PLYPD',
    'DSDPBEP86',
    'PBE0DH', 'PBEQIDH'
]


def compute_numerical_derivative_of_tensor(
        recipe, basis, derivative_repr, tensor_func, frequency=None, dry_run=False, force_choice=None, **kwargs):
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
    :param dry_run: do not fill the tensor or perform Romberg analysis
    :type dry_run: bool
    :param kwargs: args passed to tensor_func
    :type kwargs: dict
    :param force_choice: force the choice in the Romberg triangle
    :type force_choice: tuple
    :return: tensor
    :rtype: qcip_tools.derivatives.Tensor, dict
    """

    return derivatives.compute_numerical_derivative_of_tensor(
        basis,
        derivative_repr,
        tensor_func, recipe['k_max'],
        recipe['min_field'],
        recipe['ratio'],
        accuracy_level=recipe['accuracy_level'],
        frequency=frequency,
        dry_run=dry_run,
        force_choice=force_choice,
        recipe=recipe,  # kwarg
        **kwargs)


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
