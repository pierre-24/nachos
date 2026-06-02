import os
import sys
import numpy

from typing import Optional, TextIO

from qcip_tools import derivatives
from qcip_tools.chemistry_files import chemistry_datafile

from nachos.core import compute_numerical_derivative_of_tensor, fancy_output_derivative, \
    fancy_output_component_of_derivative

from nachos.core.files import Recipe


class BadBaking(Exception):
    pass


def _equal_molecules_or_raise(mol1: object, mol2: object) -> None:
    """Verify two molecular geometries have identical atomic structure.

    Args:
        mol1: First molecular geometry.
        mol2: Second molecular geometry.

    Raises:
        BadBaking: If atomic symbols differ between geometries.
    """

    if [a.symbol for a in mol1] != [a.symbol for a in mol2]:
        raise BadBaking('not the same geometries: atomic symbols are different')


class Baker:
    """Perform numerical differentiation on quantum chemistry results.

    This class computes finite-difference derivatives from computed molecular properties
    using Romberg extrapolation for improved accuracy.

    Args:
        recipe: Recipe object defining differentiation parameters.
        storage: ComputationalResults object containing computed values.
        directory: Working directory for file operations.
        original_cf: Optional existing chemistry datafile to append derivatives to.
    """

    def __init__(
            self, recipe: Recipe, storage, directory: str = '.',
            original_cf: Optional[chemistry_datafile.ChemistryDataFile] = None
    ):
        self.recipe = recipe

        if not os.path.isdir(directory):
            raise BadBaking('{} is not a directory'.format(directory))

        self.directory = directory
        self.storage = storage

        if original_cf is not None:
            _equal_molecules_or_raise(self.recipe.geometry, original_cf.molecule)
            self.original_cf = original_cf
        else:
            self.original_cf = chemistry_datafile.ChemistryDataFile.from_molecule(
                self.recipe.geometry, 'nachos ND result'
            )

        if self.storage.check() != ([], []):
            raise BadBaking('The storage (h5 file) does not fulfill the recipe!')

    def bake(
            self, only: list | None = None, out: TextIO = sys.stdout, verbosity_level: int = 0,
            copy_zero_field_basis: bool = False, force_choice: tuple | None = None
    ):
        """Compute numerical derivatives from stored quantum chemistry results.

        Performs finite-difference numerical differentiation using Romberg extrapolation
        on computed derivatives, with optional filtering and customization.

        Args:
            only: List of (derivative, level) tuples to compute (None = all).
            out: File-like object for output messages (default: sys.stdout).
            verbosity_level: Verbosity level (0=silent, 1+=verbose with increasing detail).
            copy_zero_field_basis: Copy zero-field results to all field configurations.
            force_choice: Force specific Romberg triangle choice (advanced).

        Returns:
            ChemistryDataFile object with computed derivatives appended.
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
        f = self.original_cf
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
                            frequency=freq,
                            force_choice=force_choice)
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
    def make_uncertainty_tensor(
            romberg_triangles: dict, initial_derivative: derivatives.Derivative,
            diff_derivative: derivatives.Derivative, frequency: str | float
    ) -> derivatives.Tensor:
        """Compute error estimates from Romberg extrapolation triangles.

        Creates a tensor of uncertainties by extracting convergence error estimates
        from Romberg triangles used in numerical differentiation.

        Args:
            romberg_triangles: Dictionary mapping component indices to Romberg triangles.
            initial_derivative: Starting derivative before differentiation.
            diff_derivative: Differentiation derivative to apply.
            frequency: Frequency value (for frequency-dependent properties).

        Returns:
            Tensor object containing uncertainty estimates.
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
            recipe: 'Recipe',
            initial_derivative: derivatives.Derivative,
            diff_derivative: derivatives.Derivative,
            final_result: derivatives.Tensor,
            romberg_triangles: dict,
            tensor_access,
            out: TextIO = sys.stdout,
            verbosity_level: int = 0) -> None:
        """Display detailed computation information and validation statistics.

        Outputs information about the numerical differentiation computation at various
        verbosity levels, including Romberg triangles, Kleinman conditions, and uncertainty estimates.

        Verbosity levels:
            0: Silent
            1: Output final tensor
            2: Include Romberg triangles and best value selection
            3+: Include decision process and convergence details

        Args:
            recipe: Recipe defining differentiation parameters.
            initial_derivative: Starting derivative.
            diff_derivative: Applied differentiation.
            final_result: Computed final derivative tensor.
            romberg_triangles: Romberg extrapolation triangles for each component.
            tensor_access: Function to access tensor components from storage.
            out: File-like object for output (default: sys.stdout).
            verbosity_level: Verbosity level (0-3+).
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

                        out.write('\n------------------------------------------------------\n')
                        out.write(' F          V(F)                  V(F)-V(0)\n')
                        out.write('------------------------------------------------------\n')
                        zero_field_val = tensor_access(
                            [0] * len(field), 0, initial_derivative, b_coo, final_result.frequency, recipe)

                        for i, c in enumerate(all_fields):
                            k = i - recipe['k_max']

                            field_val = 0
                            if k != 0:
                                field_val = recipe['min_field'] * recipe['ratio'] ** (abs(k) - 1) * (-1 if k < 0 else 1)

                            val = tensor_access(
                                c, 0, initial_derivative, b_coo, final_result.frequency, recipe)
                            dV = val - zero_field_val
                            out.write(
                                '{: .7f} {: .14e} {: .14e}\n'.format(
                                    field_val,
                                    val, dV))

                        out.write('------------------------------------------------------\n\n')

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

            if verbosity_level >= 2:
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


def project_geometrical_derivatives(
        recipe: 'Recipe', datafile: object, mass_weighted_hessian: object, out: TextIO = sys.stdout,
        verbosity_level: int = 0
) -> None:
    """Project geometrical derivatives onto normal modes.

    Converts geometrical derivatives from Cartesian to normal mode coordinates using
    mass-weighted Hessian transformation, useful for vibrational analysis.

    Args:
        recipe: Recipe defining differentiation parameters.
        datafile: ChemistryDataFile containing computed derivatives.
        mass_weighted_hessian: Mass-weighted Hessian for coordinate transformation.
        out: File-like object for output messages (default: sys.stdout).
        verbosity_level: Verbosity level for informational output.

    Raises:
        ValueError: If Hessian dimensions don't match recipe degrees of freedom.
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


def __project_tensor(data: object, mwh: object) -> object:
    """Project derivative tensor onto normal mode coordinates.

    Args:
        data: Derivative tensor in Cartesian coordinates.
        mwh: Mass-weighted Hessian for transformation.

    Returns:
        Derivative tensor projected onto normal modes.
    """

    return data.project_over_normal_modes(mwh)


def __output_nm_derivatives(
        recipe: 'Recipe', final_result: object, out: TextIO = sys.stdout, verbosity_level: int = 0,
        trans_plus_rot_dof: int = 0
) -> None:
    """Display normal mode projected derivatives with formatting.

    Args:
        recipe: Recipe defining parameters.
        final_result: Projected derivative tensor to display.
        out: File-like object for output (default: sys.stdout).
        verbosity_level: Verbosity level for output.
        trans_plus_rot_dof: Number of translational/rotational degrees of freedom to skip.
    """
    if verbosity_level >= 1:
        out.write('\n*** projected ')
        out.write(fancy_output_derivative(final_result.representation, final_result.frequency))
        out.write('\n')
        out.write(final_result.to_string(
            molecule=recipe.geometry, threshold=1e-8, skip_trans_plus_rot_dof=trans_plus_rot_dof))
