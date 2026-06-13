import os
import glob
import math
from typing import TextIO

import numpy
import sys

from qcip_tools import quantities, derivatives, derivatives_e
from qcip_tools.molecule import Molecule
from qcip_tools.chemistry_files import helpers, PropertyNotPresent, PropertyNotDefined

from nachos.core import files, preparing, GAUSSIAN_DOUBLE_HYBRIDS
from nachos.qcip_tools_ext import gaussian, qchem, dalton  # noqa
from nachos.core.files import ComputationalResults


class BadCooking(Exception):
    pass


class Cooker:
    """Extract computed properties from quantum chemistry calculation results.

    This class processes output files from quantum chemistry packages and collects
    computed derivatives for storage in a computational results file.

    Args:
        recipe: Recipe object defining the differentiation parameters.
        directory: Working directory where storage will be written.
    """

    def __init__(self, recipe, directory='.'):
        self.recipe = recipe

        if not os.path.isdir(directory):
            raise BadCooking('{} is not a directory'.format(directory))

        self.directory = directory
        self.fields_needed_by_recipe = preparing.fields_needed_by_recipe(self.recipe)
        self.fields_needed = [a[0] for a in self.fields_needed_by_recipe]

    def cook(
            self, directories: list[str], out: TextIO = sys.stdout, verbosity_level: int = 0,
            use_gaussian_logs: bool = False
    ) -> 'ComputationalResults':
        """Process quantum chemistry output files and collect computed derivatives.

        Searches directories for calculation output files, extracts computed properties
        (energies, derivatives), and stores them in a ComputationalResults object.

        Args:
            directories: List of directories to search for QM results.
            out: File-like object for output messages (default: sys.stdout).
            verbosity_level: Verbosity level (0=silent, 1=verbose).
            use_gaussian_logs: Use Gaussian LOG files instead of FCHK (not recommended).

        Returns:
            ComputationalResults object containing all extracted derivatives.
        """

        storage = files.ComputationalResults(self.recipe, directory=self.directory)

        if self.recipe['flavor'] == 'gaussian':
            look_for = ['*.fchk']
            if use_gaussian_logs:
                look_for = ['*.log']
        elif self.recipe['flavor'] == 'dalton':
            look_for = ['*.tar.gz']
            if any(a[0] == 'G' for a in self.recipe.bases()):
                look_for.append('*.out')
        elif self.recipe['flavor'] == 'qchem':
            look_for = ['*.log']
        else:
            look_for = []

        for directory in directories:
            if not os.path.isdir(directory):
                raise BadCooking('{} is no directory!'.format(directory))

            for line in look_for:
                for i in glob.glob('{}/{}'.format(directory, line)):
                    if verbosity_level >= 1:
                        out.write('* cooking with {} ... '.format(i))

                    with open(i) as f:
                        try:
                            fx = helpers.open_chemistry_file(f)
                            obtained = self.cook_from_file(fx, i, storage)
                            if verbosity_level >= 1:
                                out.write(
                                    '({}) ... ok\n'.format(', '.join(obtained))
                                    if len(obtained) != 0 else 'empty\n')

                        except helpers.ProbablyNotAChemistryFile:
                            if verbosity_level >= 1:
                                out.write('skipped\n')
                            continue

        return storage

    def cook_from_file(self, f: object, name: str, storage: 'ComputationalResults') -> list[str]:
        """Extract derivatives from a single quantum chemistry output file.

        Processes a chemistry file to extract computed energies, electrical derivatives,
        and geometrical derivatives, storing results in the provided storage object.

        Args:
            f: Chemistry file object opened for reading.
            name: Path to the file being processed.
            storage: ComputationalResults object to store extracted data.

        Returns:
            List of derivative identifiers that were successfully extracted.
        """

        def almost_the_same(a, b, threshold=1e-3):
            return math.fabs(a - b) < threshold

        obtained = []

        # catch field
        if self.recipe['type'] == 'F':
            try:
                real_fields = f.property('n:input_electric_field')[1:4]

                if f.file_type in ['GAUSSIAN_FCHK', 'GAUSSIAN_LOG']:
                    # Gaussian computes actually for the opposite field!
                    real_fields = list(-x for x in real_fields)

            except PropertyNotPresent:
                raise BadCooking('F derivative but not electric field for {}'.format(name))
        else:
            try:
                real_fields = Cooker.real_fields_from_geometry(
                    self.recipe.geometry, f.property('molecule'), threshold=.5 * self.recipe['min_field'])
            except (PropertyNotDefined, PropertyNotPresent) as e:
                raise BadCooking('G derivative but unable to get geometry from {} ({})'.format(name, str(e)))

        fields = Cooker.real_fields_to_fields(real_fields, self.recipe['min_field'], self.recipe['ratio'])

        if fields not in self.fields_needed:
            return []  # not part of this calculation

        level = next(a[1] for a in self.fields_needed_by_recipe if a[0] == fields)
        derivatives_per_level = self.recipe.bases(level_min=level)
        derivatives_in_level = [a[0] for a in derivatives_per_level]

        if derivatives.Derivative() in derivatives_in_level and f.file_type != 'DALTON_LOG':
            try:
                energies = f.property('computed_energies')

                if f.file_type in ['GAUSSIAN_FCHK', 'QCHEM_LOG', 'GAUSSIAN_LOG']:
                    if self.recipe['method'] == 'SCS-MP2':
                        # use the parameters from S. Grimme. J. Chem. Phys. 118, 9095 (2003).
                        energy = energies['HF']
                        sc_energies = f.property('n:spin_components_e2')
                        energy += 1 / 3 * (sc_energies['aa'] + sc_energies['bb'])  # p_T * E_T
                        energy += 6 / 5 * sc_energies['ab']  # p_S * E_S
                    else:
                        key = self.recipe['method']
                        if key == 'DFT':
                            key = 'SCF/DFT'
                            if self.recipe['flavor_extra']['XC'] in GAUSSIAN_DOUBLE_HYBRIDS:
                                key = 'MP2'  # Gaussian puts the result in the MP2 field instead of SCF/DFT.
                        energy = energies[key]  # tries to catch the energy for the correct method
                    obtained.append('energy:' + self.recipe['method'])
                else:  # ?!?
                    energy = energies['total']
                    obtained.append('energy')

                storage.add_result(
                    fields,
                    '',
                    derivatives.Tensor('', components=numpy.array((energy,))),
                    allow_replace=True)

            except (PropertyNotPresent, PropertyNotDefined, KeyError):
                pass

        try:
            electrical_derivatives = f.property('electrical_derivatives')
            for d in electrical_derivatives:
                if d in derivatives_in_level:
                    e_deriv = {}
                    for a, b in electrical_derivatives[d].items():
                        freq = a
                        if a != 'static':
                            for freq_ir in self.recipe['frequencies']:
                                if almost_the_same(
                                        derivatives_e.convert_frequency_from_string(freq_ir),
                                        derivatives_e.convert_frequency_from_string(freq)):
                                    freq = freq_ir
                                    break

                        e_deriv[freq] = b

                    obtained.append(d)

                    storage.add_result(
                        fields,
                        d,
                        e_deriv,
                        # NOTE: frequency calculations compute electrical deriv as well, so replace them if any:
                        allow_replace=d in ['F', 'FF', 'FFF'])
        except (PropertyNotPresent, PropertyNotDefined):
            pass

        try:
            geometrical_derivatives = f.property('geometrical_derivatives')
            for d in geometrical_derivatives:
                if d in derivatives_in_level:
                    storage.add_result(fields, d, geometrical_derivatives[d])
        except (PropertyNotPresent, PropertyNotDefined):
            pass
        return obtained

    @staticmethod
    def real_fields_to_fields(real_field: list, min_field: float, ratio: float) -> list[int]:
        """Convert real field values to discretized field indices.

        Transforms continuous field values into discrete indices based on geometric
        progression with specified minimum field and ratio.

        Args:
            real_field: List of real-valued field strengths.
            min_field: Minimum field value (first non-zero index).
            ratio: Geometric progression ratio for field discretization.

        Returns:
            List of integer field indices (0 for zero values).
        """
        fields = []
        for val in real_field:
            if val == .0:
                fields.append(0)
            else:
                sgn = -1 if val < 0 else 1
                i = round(math.log(math.fabs(val) / min_field) / math.log(ratio)) + 1
                fields.append(sgn * i)

        return fields

    @staticmethod
    def real_fields_from_geometry(
            geometry: Molecule, deformed_geometry: Molecule, threshold: float = 1e-4) -> list[float]:
        """Calculate real-valued field displacements from geometry deformation.

        Compares two molecular geometries to extract field-dependent atomic displacements,
        useful for identifying field values from geometrically distorted structures.

        Args:
            geometry: Reference molecular geometry.
            deformed_geometry: Deformed molecular geometry.
            threshold: Minimum displacement magnitude to consider non-zero (default: 1e-4).

        Returns:
            List of real-valued field displacements in atomic units.

        Raises:
            ValueError: If geometries have different atoms or structures.
        """
        if len(geometry) != len(deformed_geometry) or str(geometry) != str(deformed_geometry):
            raise ValueError('geometries does not contain the same number of atoms')

        real_fields = [0] * (3 * len(geometry))

        for index, a in enumerate(deformed_geometry):
            for i in range(3):
                diff_value = (a.position[i] - geometry[index].position[i]) / quantities.AuToAngstrom
                real_fields[index * 3 + i] = 0 if math.fabs(diff_value) < threshold else diff_value

        return real_fields
