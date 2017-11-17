import os
import sys

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

        for initial_derivative, level in bases:
            for diff_order in range(1, level + 1):
                diff_derivative = derivatives.Derivative(self.recipe['type'] * diff_order, spacial_dof=dof)
                final_derivative = initial_derivative.differentiate(diff_derivative.representation())

                if derivatives.is_electrical(initial_derivative) or self.recipe['type'] == 'F':
                    results_per_frequency = {}
                    if 'D' in initial_derivative.representation():
                        for freq in self.recipe['frequencies']:
                            r, trigs = compute_numerical_derivative_of_tensor(
                                self.recipe,
                                initial_derivative,
                                diff_derivative,
                                self.storage.tensor_element_access,
                                frequency=freq)
                            results_per_frequency[freq] = r.components

                            Baker.output_information(self.recipe, initial_derivative, r, trigs, out, verbosity_level)
                    else:
                        r, trigs = compute_numerical_derivative_of_tensor(
                            self.recipe,
                            initial_derivative,
                            diff_derivative,
                            self.storage.tensor_element_access,
                            frequency='static')
                        results_per_frequency['static'] = r.components
                        Baker.output_information(self.recipe, initial_derivative, r, trigs, out, verbosity_level)

                    f.derivatives[final_derivative.representation()] = results_per_frequency
                else:
                    r, trigs = compute_numerical_derivative_of_tensor(
                        self.recipe, initial_derivative, diff_derivative, self.storage.tensor_element_access)
                    f.derivatives[final_derivative.representation()] = r.components

                    Baker.output_information(self.recipe, initial_derivative, r, trigs, out, verbosity_level)
        return f

    @staticmethod
    def output_information(
            recipe, initial_derivative, final_result, romberg_triangles, out=sys.stdout, verbosity_level=0):
        """Output information about what was computed

        .. note::

            If verbosity level is:

            - **<=0:** nothing happen ;
            - **>0:** output final tensor ;
            - **>1:** output Romberg triangle and best value ;
            - **>2:** output decision process to find best value in Romberg triangle.

            Therefore, it triggers once again the computation of the best value in Romberg triangle
            if verbosity_level is > 1.

        :param recipe: the corresponding recipe
        :type recipe: nachos.core.files.Recipe
        :param out: output
        :type out: file
        :param initial_derivative: starting point
        :type initial_derivative: qcip_tools.derivatives.Derivative
        :param final_result: what was computed
        :type final_result: qcip_tools.derivatives.Tensor
        :param romberg_triangles: the different Romberg triangles
        :type romberg_triangles: collections.OrderedDict
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

                    for b_coo in romberg_triangles[d_coo]:
                        if initial_derivative != '':
                            out.write('\n# component {} of {}:\n'.format(
                                fancy_output_component_of_derivative(initial_derivative, b_coo, recipe.geometry),
                                basis_name))

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
