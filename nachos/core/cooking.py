import collections
import copy
import math
import os

import numpy
from qcip_tools import numerical_differentiation, derivatives, quantities, derivatives_e
from qcip_tools.chemistry_files import gaussian


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
        if recipe['type'] == 'electric':
            # because E(-F) = E0 - µ.F + 1/2*a.F² + ..
            #         µ(-F) = µ0 - a.F + 1/2*b.F² + ...
            # (and why not G and GG ?)
            if b_repr == '':
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
                            tensor.components[inv_basis_coo][inv_derivative_coo] = val[1]
                        else:
                            tensor.components[inv_derivative_coo] = val[1]

    if not dry_run:
        return tensor, romberg_triangles


def fields_needed_by_recipe(recipe):
    """

    :param recipe: recipe
    :type recipe: nachos.core.files.Recipe
    :rtype: list
    """

    diffs = recipe.maximum_derivatives()

    fields_needed = [
        [0] * (3 if recipe['type'] == 'F' else recipe.dof)  # with zero field
    ]

    fields_needed_with_level = [
        (fields_needed[0], 1)
    ]

    def collect_fields(fields, *args, **kwargs):
        f = list(int(a) for a in fields)
        if f not in fields_needed:
            fields_needed.append(f)
            fields_needed_with_level.append((f, kwargs['level']))
        return .0

    e = derivatives.Derivative('')

    for d in diffs:
        compute_numerical_derivative_of_tensor(
            recipe, e, d, collect_fields, dry_run=True, level=d.order())

    return fields_needed_with_level


class BadCooking(Exception):
    pass


class Cooker:
    """Cooking class to generate the input files

    :param recipe: a recipe
    :type recipe: nachos.core.files.Recipe
    """

    def __init__(self, recipe, directory='.'):
        self.recipe = recipe

        if not os.path.isdir(directory):
            raise BadCooking('{} is not a directory'.format(directory))

        self.directory = directory
        self.fields_needed = fields_needed_by_recipe(self.recipe)

    def cook(self):
        """Create the different input files in the directory"""

        return getattr(self, 'cook_{}_inputs'.format(self.recipe['flavor']))()

    def cook_gaussian_inputs(self):
        """Create inputs for gaussian
        """

        base_m = False
        counter = 0
        files_created = 0

        # try to open custom basis set if any
        gen_basis_set = []
        if self.recipe['basis_set'] == 'gen':
            path = self.recipe['flavor_extra']['gen_basis']
            if not os.path.exists(path):
                raise BadCooking('gen basis file {} cannot be opened'.format(path))

            gbs = gaussian.BasisSet()
            try:
                with open(path) as f:
                    gbs.read(f)
            except gaussian.BasisSetFormatError as e:
                raise BadCooking('error while opening custom basis set ({}): {}'.format(path, str(e)))

            try:
                gen_basis_set = gbs.to_string(for_molecule=self.recipe.geometry).splitlines()[1:]
            except Exception as e:
                raise BadCooking('error while using custom basis set ({}) : {}'.format(path, e))

        for fields, level in self.fields_needed:
            counter += 1

            compute_polar = False
            compute_polar_with_freq = False
            compute_G = False
            compute_GG = False

            bases = self.recipe.bases(level_min=level)

            for b, l in bases:
                basis = b.representation()
                if not compute_G and basis == 'G':
                    compute_G = True
                if not compute_GG and basis == 'GG':
                    compute_GG = True
                if not compute_polar and 'FF' in basis:
                    compute_polar = True
                if not compute_polar_with_freq and 'D' in basis:
                    compute_polar = True
                    compute_polar_with_freq = True

            compute_polar_and_G = compute_polar and (compute_G or compute_GG)

            fi = gaussian.Input()
            real_fields = Cooker.real_fields(fields, self.recipe['min_field'], self.recipe['ratio'])

            if self.recipe['type'] == 'G':
                fi.molecule = Cooker.deform_geometry(self.recipe.geometry, real_fields)
            else:
                fi.molecule = self.recipe.geometry

            if not base_m and fields == [0] * len(fields):
                fi.title = 'base'
                base_m = True
            else:
                fi.title = 'field({})='.format(level) + \
                    ','.join(Cooker.nonzero_fields(fields, self.recipe.geometry, self.recipe['type']))

            fi.options['nprocshared'] = self.recipe['flavor_extra']['procs']
            fi.options['mem'] = self.recipe['flavor_extra']['memory']
            fi.options['chk'] = 'xxx'

            # input card
            input_card = [
                '#P {} {} nosym{}'.format(
                    self.recipe['method'] if self.recipe['method'] != 'DFT' else self.recipe['flavor_extra']['XC'],
                    self.recipe['basis_set'],
                    ' field=read' if self.recipe['type'] == 'F' else ''),
                'scf=(Conver={c},NoVarAcc,MaxCyc={m},vshift={v}) IOP(9/6={m},9/9={cc})'.format_map(
                    {
                        'c': self.recipe['flavor_extra']['convergence'],
                        'cc': self.recipe['flavor_extra']['cc_convergence'],
                        'm': self.recipe['flavor_extra']['max_cycles'],
                        'v': self.recipe['flavor_extra']['vshift']
                    })
            ]

            if self.recipe['flavor_extra']['extra_keywords']:
                input_card.extend(self.recipe['flavor_extra']['extra_keywords'])

            fi.input_card = input_card

            # other blocks
            if self.recipe['basis_set'] == 'gen':
                fi.other_blocks.append(gen_basis_set)

            if self.recipe['type'] == 'F':
                fi.other_blocks.append(['\t'.join(['{: .10f}'.format(a) for a in real_fields])])

            if self.recipe['flavor_extra']['extra_sections']:
                fi.other_blocks.extend(self.recipe['flavor_extra']['extra_sections'])

            # write files
            if compute_polar:
                extra_line = 'polar{} cphf=(conver={}{})'.format(
                    '=dcshg' if 'FDD' in bases else '',
                    self.recipe['flavor_extra']['cphf_convergence'],
                    ',rdfreq' if compute_polar_with_freq else ''
                )
                fi.input_card.append(extra_line)

                if compute_polar_with_freq:
                    fi.other_blocks.insert(0, [
                        '{}'.format(derivatives_e.convert_frequency_from_string(a)) for a in
                        self.recipe['frequencies']])
                with open('{}/{}_{:04d}{}.com'.format(
                        self.directory, self.recipe['name'], counter, 'a' if compute_polar_and_G else ''), 'w') as f:
                    fi.write(f)
                    files_created += 1
                fi.input_card.pop(-1)
                fi.other_blocks.pop(0)

            if compute_polar_and_G or not compute_polar:
                if compute_GG:
                    fi.input_card.append('freq')
                elif compute_G:
                    fi.input_card.append('force')

                with open('{}/{}_{:04d}{}.com'.format(
                        self.directory, self.recipe['name'], counter, 'b' if compute_polar_and_G else ''), 'w') as f:
                    fi.write(f)
                    files_created += 1

        return files_created

    def cook_dalton_inputs(self):
        """Create inputs for dalton
        """

        for fields in self.fields_needed:
            real_fields = Cooker.real_fields(fields, self.recipe['min_field'], self.recipe['ratio'])
            print(fields, real_fields)

    @staticmethod
    def deform_geometry(geometry, fields, geometry_in_angstrom=True):
        """Create an input for gaussian

        :param fields: Differentiation field
        :type fields: list
        :param geometry: geometry do deform
        :type geometry: qcip_tools.molecule.Molecule
        :param geometry_in_angstrom: indicate wheter the geometry is given in Angstrom or not
            (because the field is obviously given in atomic units)
        :type geometry_in_angstrom: bool
        :rtype: qcip_tools.molecule.Molecule
        """

        deformed = copy.deepcopy(geometry)

        for index, atom in enumerate(deformed):
            atom.position += numpy.array(
                [fields[index * 3 + i] * (quantities.AuToAngstrom if geometry_in_angstrom else 1.) for i in range(3)]
            )

        return deformed

    @staticmethod
    def real_fields(fields, min_field, ratio):
        """Return the "real value" of the field applied

        :param fields: input field (in term of ak)
        :type fields: list
        :param min_field: minimal field
        :type min_field: float
        :param ratio: ratio
        :type ratio: float
        :rtype: list
        """

        return [min_field * numerical_differentiation.ak_shifted(ratio, _) for _ in fields]

    @staticmethod
    def nonzero_fields(fields, geometry, t):
        return [
            '{}({:+g}{})'.format(
                '{}{}'.format(
                    geometry[int(math.floor(i / 3))].symbol, int(math.floor(i / 3) + 1))
                if t == 'G' else 'F',
                e,
                derivatives.COORDINATES[i % 3]
            ) for i, e in enumerate(fields) if e != 0]


class Baker:
    """Baking class to retrieve the information out of the calculation results

    :param recipe: a recipe
    :type recipe: nachos.core.files.Recipe"""

    def __init__(self, recipe):
        self.recipe = recipe
