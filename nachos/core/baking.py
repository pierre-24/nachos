import os
import sys
import numpy

from qcip_tools import derivatives
from qcip_tools.chemistry_files import chemistry_datafile

from nachos.core import compute_numerical_derivative_of_tensor, fancy_output_derivative, \
    fancy_output_component_of_derivative


class BadBaking(Exception):
    pass


class Baker:
    """Baker class to finally perform the numerical differentiation

    :param recipe: a recipe
    :type recipe: nachos.core.files.Recipe
    :param storage: storage of results
    :type storage: nachos.core.files.ComputationalResults
    :param directory: working directory
    :type directory: str
    """

    def __init__(self, recipe, storage, directory='.'):
        self.recipe = recipe

        if not os.path.isdir(directory):
            raise BadBaking('{} is not a directory'.format(directory))

        self.directory = directory
        self.storage = storage

        if self.storage.check() != ([], []):
            raise BadBaking('The storage (h5 file) does not fulfill the recipe!')

    def bake(self, only=None, out=sys.stdout, verbosity_level=0, copy_zero_field_basis=False):
        """Perform the numerical differentiation

        :param only: list of derivatives to perform (None = all of them)
        :type only: list
        :param out: output if information is needed to be outputed
        :type out: file
        :param verbosity_level: how far should we print information
        :type verbosity_level: int
        :param copy_zero_field_basis: copy the basis found in zero field
        :type copy_zero_field_basis: bool
        :rtype: qcip_tools.chemistry_files.chemistry_datafile.ChemistryDataFile
        """

        if not only:
            bases = [a for a in self.recipe.bases()]
        else:
            bases = []
            for base, level in self.recipe.bases():
                try:
                    _, n_level = next(a for a in only if a[0] == base)
                    if n_level > level:
                        raise BadBaking('level of differentiation for {} ({}) is larger than available ({})'.format(
                            base if base != '' else 'energy', n_level, level))
                    if n_level < 1:
                        n_level = level

                    bases.append((base, n_level))
                except StopIteration:
                    continue

        if len(bases) == 0:
            raise BadBaking('no differentiation requested!')

        bases.sort(key=lambda a: a[1], reverse=True)
        f = chemistry_datafile.ChemistryDataFile.from_molecule(self.recipe.geometry, 'nachos ND result')
        dof = 3 * len(self.recipe.geometry)

        if copy_zero_field_basis:
            zero_field = tuple([0] * (dof if self.recipe['type'] == 'G' else 3))
            for b in self.storage.results[zero_field]:
                f.derivatives[b] = self.storage.results[zero_field][b]

        # print fields (by request of Benoit)
        if verbosity_level > 1:
            fields = []
            for i in range(0, self.recipe['k_max']):
                fields.append(self.recipe['min_field'] * self.recipe['ratio'] ** i)

            out.write('! Note:\n')
            out.write('! Type of differentiation is: {}.\n'.format(
                'electrical' if self.recipe['type'] == 'F' else 'geometrical'))
            out.write('! Minimum field is: {} (a.u), ratio is: {}, k_max is: {}.\n'.format(
                self.recipe['min_field'], self.recipe['ratio'], self.recipe['k_max']))
            out.write('! Thus, fields used (a.u.) during differentiation are: {}.\n\n'.format(
                ', '.join('{}'.format(i) for i in fields)))

        for initial_derivative, level in bases:
            for diff_order in range(1, level + 1):
                diff_derivative = derivatives.Derivative(self.recipe['type'] * diff_order, spacial_dof=dof)
                final_derivative = initial_derivative.differentiate(diff_derivative.representation())

                if derivatives.is_electrical(initial_derivative) or self.recipe['type'] == 'F':
                    results_per_frequency = {}

                    freqs = []
                    if 'D' in initial_derivative.representation():
                        freqs.extend(self.recipe['frequencies'])
                    else:
                        freqs = ['static']

                    for freq in freqs:
                        r, trigs = compute_numerical_derivative_of_tensor(
                            self.recipe,
                            initial_derivative,
                            diff_derivative,
                            self.storage.tensor_element_access,
                            frequency=freq)
                        results_per_frequency[freq] = r
                        Baker.output_information(
                            self.recipe,
                            initial_derivative,
                            diff_derivative,
                            r,
                            trigs,
                            self.storage.tensor_element_access,
                            out,
                            verbosity_level)

                    f.derivatives[final_derivative.representation()] = results_per_frequency
                else:
                    r, trigs = compute_numerical_derivative_of_tensor(
                        self.recipe,
                        initial_derivative,
                        diff_derivative,
                        self.storage.tensor_element_access)
                    f.derivatives[final_derivative.representation()] = r

                    Baker.output_information(
                        self.recipe,
                        initial_derivative,
                        diff_derivative,
                        r,
                        trigs,
                        self.storage.tensor_element_access,
                        out,
                        verbosity_level)
        return f

    @staticmethod
    def make_uncertainty_tensor(romberg_triangles, initial_derivative, diff_derivative, frequency):
        """

        :param romberg_triangles: the different Romberg triangles
        :type romberg_triangles: collections.OrderedDict
        :param initial_derivative: starting point
        :type initial_derivative: qcip_tools.derivatives.Derivative
        :param diff_derivative: differentialtion
        :type diff_derivative: qcip_tools.derivatives.Derivative
        :param frequency: the frequency
        :type frequency: str|float
        :rtype: qcip_tools.derivatives.Tensor
        """

        final_derivative = initial_derivative.differentiate(diff_derivative.representation())
        t = derivatives.Tensor(final_derivative, spacial_dof=final_derivative.spacial_dof, frequency=frequency)

        for d_coo in romberg_triangles:
            for b_coo in romberg_triangles[d_coo]:
                total_coo = list(d_coo if type(d_coo) is tuple else (d_coo, ))
                if final_derivative.basis != '':
                    total_coo.extend(b_coo if type(b_coo) is tuple else (b_coo, ))

                uncertainty = romberg_triangles[d_coo][b_coo]()[-1]

                for e in initial_derivative.inverse_smart_iterator(b_coo):
                    for ex in diff_derivative.inverse_smart_iterator(d_coo):
                        if initial_derivative.representation() != '':
                            if 'G' in diff_derivative.representation():
                                t.components[ex][e] = uncertainty
                            else:
                                t.components[e][ex] = uncertainty
                        else:
                            t.components[ex] = uncertainty

        return t

    @staticmethod
    def output_information(
            recipe,
            initial_derivative,
            diff_derivative,
            final_result,
            romberg_triangles,
            tensor_access,
            out=sys.stdout,
            verbosity_level=0):
        """Output information about what was computed

        .. note::

            If verbosity level is:

            - **<=0:** nothing happen ;
            - **>0:** output final tensor ;
            - **>1:** output Romberg triangle and best value ;
            - **>2:** output decision process to find best value in Romberg triangle.

            Therefore, it triggers once again the computation of the best value in Romberg triangle
            if verbosity_level is > 2.

        :param recipe: the corresponding recipe
        :type recipe: nachos.core.files.Recipe
        :param initial_derivative: starting point
        :type initial_derivative: qcip_tools.derivatives.Derivative
        :param final_result: what was computed
        :type final_result: qcip_tools.derivatives.Tensor
        :param romberg_triangles: the different Romberg triangles
        :type romberg_triangles: collections.OrderedDict
        :param out: output
        :type out: file
        :param verbosity_level: how far should we print information
        :type verbosity_level: int
        """

        if verbosity_level >= 1:
            basis_name = fancy_output_derivative(initial_derivative)
            out.write('*** {} derivative of {} to get {}:\n'.format(
                'geometrical' if recipe['type'] == 'G' else 'electrical',
                basis_name,
                fancy_output_derivative(final_result.representation, final_result.frequency)))

            if verbosity_level >= 2:
                out.write('** Romberg triangles:\n')
                for d_coo in romberg_triangles:
                    out.write('\n* computing {} / {}:\n'.format(
                        fancy_output_derivative(initial_derivative, final_result.frequency),
                        ' '.join(
                            derivatives.representation_to_operator(recipe['type'], a, recipe.geometry) for a in d_coo)
                    ))

                    # generate fields
                    if 'G' in diff_derivative.representation():
                        field = [0] * diff_derivative.spacial_dof
                    else:
                        field = [0] * 3

                    for b in d_coo:
                        field[b] = 1

                    all_fields = [[0] * len(field)]
                    for i in range(1, recipe['k_max'] + 1):
                        all_fields.append(list(x * i for x in field))
                        all_fields.insert(0, list(-x * i for x in field))

                    for b_coo in romberg_triangles[d_coo]:
                        if initial_derivative != '':
                            out.write('\n# component {} of {}:\n'.format(
                                fancy_output_component_of_derivative(initial_derivative, b_coo, recipe.geometry),
                                basis_name))

                        out.write('\n----------------------------------------------\n')
                        out.write(' F          V(F)              V(F)-V(0)\n')
                        out.write('----------------------------------------------\n')
                        zero_field_val = tensor_access(
                            [0] * len(field), 0, initial_derivative, False, b_coo, final_result.frequency, recipe)

                        for i, c in enumerate(all_fields):
                            k = i - recipe['k_max']

                            field_val = 0
                            if k != 0:
                                field_val = recipe['min_field'] * recipe['ratio'] ** (abs(k) - 1) * (-1 if k < 0 else 1)

                            val = tensor_access(
                                c, 0, initial_derivative, False, b_coo, final_result.frequency, recipe)
                            dV = val - zero_field_val
                            out.write(
                                '{: .7f} {: .10e} {: .10e}\n'.format(
                                    field_val,
                                    val, dV))

                        out.write('----------------------------------------------\n\n')

                        romberg_triangle = romberg_triangles[d_coo][b_coo]
                        out.write(romberg_triangle.romberg_triangle_repr(with_decoration=True))

                        vals = romberg_triangle.find_best_value(verbose=verbosity_level >= 3, out=out)
                        if verbosity_level >= 3:
                            out.write('\n')

                        out.write('({}) = {:.5e} (error = {:.5e})\n'.format(
                            ','.join(str(a) for a in vals[0]), vals[1], vals[2]))

            if verbosity_level >= 2:
                out.write('\n** Final result:\n')

            out.write(final_result.to_string(molecule=recipe.geometry))
            out.write('\n')

            if verbosity_level >= 2 and initial_derivative != '':
                out.write('\n** Checking Kleinman conditions:\n')
                for i in final_result.representation.smart_iterator():
                    values = list(
                        final_result.components[j] for j in final_result.representation.inverse_smart_iterator(i))
                    if len(values) > 1:
                        out.write('- {}: '.format(fancy_output_component_of_derivative(final_result.representation, i)))
                        out.write('{: .5e} ± {:.5e}'.format(numpy.mean(values), numpy.std(values)))
                        out.write('\n')

                out.write('\n** Estimated uncertainties:\n')
                out.write('*** Values:\n')

                u = Baker.make_uncertainty_tensor(
                    romberg_triangles, initial_derivative, diff_derivative, final_result.frequency)

                out.write(u.to_string(molecule=recipe.geometry, threshold=1e-8))
                out.write('\n')

                out.write('*** Ratio (%):\n')
                ru = derivatives.Tensor(
                    u.representation, components=(u.components / final_result.components) * 100,
                    spacial_dof=u.spacial_dof,
                    frequency=u.frequency
                )

                out.write(ru.to_string(molecule=recipe.geometry))
                out.write('\n')


def project_geometrical_derivatives(recipe, datafile, mass_weighted_hessian, out=sys.stdout, verbosity_level=0):
    """Project geometrical derivatives, if any

    :param recipe: the recipe
    :type recipe: nachos.core.files.Recipe
    :param datafile: the data file with derivatives
    :type datafile: qcip_tools.chemistry_files.chemistry_datafile.ChemistryDataFile
    :param mass_weighted_hessian: the mass weighted hessian
    :type mass_weighted_hessian: qcip_tools.derivatives_g.MassWeightedHessian
    :param out: output
    :type out: file
    :param verbosity_level: how far should we print information
    :type verbosity_level: int
    """

    if mass_weighted_hessian.dof != recipe.dof:
        raise ValueError('displacement shape does not match')

    for basis, level in recipe.bases():
        b_repr = basis.representation()

        if derivatives.is_geometrical(b_repr) or recipe['type'] == 'G':
            for lvl in range(0, level + 1):
                if lvl == 0:
                    derivative = basis
                else:
                    derivative = basis.differentiate(recipe['type'] * lvl)
                d_repr = derivative.representation()
                n_repr = d_repr.replace('G', 'N')
                if d_repr in datafile.derivatives and n_repr not in datafile.derivatives:
                    if derivatives.is_electrical(derivative):
                        x = {}
                        for freq in datafile.derivatives[d_repr]:
                            r = __project_tensor(
                                datafile.derivatives[d_repr][freq], mass_weighted_hessian)
                            __output_nm_derivatives(recipe, r, out, verbosity_level, datafile.trans_plus_rot_dof)
                            x[freq] = r
                    else:
                        r = __project_tensor(datafile.derivatives[d_repr], mass_weighted_hessian)
                        __output_nm_derivatives(recipe, r, out, verbosity_level, datafile.trans_plus_rot_dof)
                        x = r

                    datafile.derivatives[n_repr] = x


def __project_tensor(data, mwh):
    """Project tensor over displacements to get their normal derivative equivalent

    :param data: the data
    :type data: qcip_tools.derivatives.Tensor
    :param mwh: mass weighted hessian
    :type mwh: qcip_tools.derivatives_g.MassWeightedHessian
    :rtype: qcip_tools.derivatives.Tensor
    """

    return data.project_over_normal_modes(mwh)


def __output_nm_derivatives(recipe, final_result, out, verbosity_level=0, trans_plus_rot_dof=0):
    if verbosity_level >= 1:
        out.write('\n*** projected ')
        out.write(fancy_output_derivative(final_result.representation, final_result.frequency))
        out.write('\n')
        out.write(final_result.to_string(
            molecule=recipe.geometry, threshold=1e-8, skip_trans_plus_rot_dof=trans_plus_rot_dof))
