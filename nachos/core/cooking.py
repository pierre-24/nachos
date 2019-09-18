import os
import glob
import math
import numpy
import sys

from qcip_tools import quantities, derivatives, derivatives_g, derivatives_e
from qcip_tools.chemistry_files import helpers, gaussian, PropertyNotPresent, PropertyNotDefined, dalton

from nachos.core import files, preparing


class BadCooking(Exception):
    pass


@gaussian.FCHK.define_property('n:input_electric_field')
def gaussian__FCHK__get_input_electric_field(obj, *args, **kwargs):
    """Get the input electric field

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.FCHK
    :rtype: dict
    """

    if 'External E-field' in obj:
        field = obj.get('External E-field')
        return field
    else:
        raise PropertyNotPresent('n:input_electric_field')


@dalton.Output.define_property('geometrical_derivatives')
def dalton__output__get_gradient(obj, *args, **kwargs):
    """Get the cartesian gradient out of a dalton calculation output

    (redefined so that it only fetch the gradient)

    :param obj: object
    :type obj: qcip_tools.chemistry_files.dalton.Output
    :rtype: dict
    """

    geometrical_derivatives = {}

    dof = 3 * len(obj.molecule)
    trans_plus_rot = 5 if obj.molecule.linear() else 6

    # CC
    x = obj.search('Molecular gradient', into='CC')
    if x != -1:
        gradient = derivatives_g.BaseGeometricalDerivativeTensor(
            representation='G', spacial_dof=dof, trans_plus_rot=trans_plus_rot)

        for i in range(len(obj.molecule)):
            line = obj.lines[x + 3 + i].split()
            gradient.components[i * 3:(i + 1) * 3] = [float(a) for a in line[1:]]
        geometrical_derivatives['G'] = gradient

    # general !
    if obj.chunk_exists('ABACUS'):
        x = obj.search('Molecular gradient', into='ABACUS')
        if x != -1:
            gradient = derivatives_g.BaseGeometricalDerivativeTensor(
                representation='G', spacial_dof=dof, trans_plus_rot=trans_plus_rot)

            for i in range(len(obj.molecule)):
                line = obj.lines[x + 3 + i].split()
                gradient.components[i * 3:(i + 1) * 3] = [float(a) for a in line[1:]]
            geometrical_derivatives['G'] = gradient

    if not geometrical_derivatives:
        raise PropertyNotPresent('geometrical_derivatives')

    return geometrical_derivatives


class Cooker:
    """Cooker class to retrieve the information out of the calculation results

    :param recipe: a recipe
    :type recipe: nachos.core.files.Recipe
    :param directory: working directory, where storage will be written
    :type directory: str
    """

    def __init__(self, recipe, directory='.'):
        self.recipe = recipe

        if not os.path.isdir(directory):
            raise BadCooking('{} is not a directory'.format(directory))

        self.directory = directory
        self.fields_needed_by_recipe = preparing.fields_needed_by_recipe(self.recipe)
        self.fields_needed = [a[0] for a in self.fields_needed_by_recipe]

    def cook(self, directories, out=sys.stdout, verbosity_level=0):
        """Cook files in directories, all together in a storage file

        :param directories: directories where QM results should be looked for
        :type directories: list of str
        :param out: outpout of eventual information
        :type out: file
        :param verbosity_level: wheter to write information or not
        :type verbosity_level: bool
        :rtype: nachos.core.files.ComputationalResults
        """

        storage = files.ComputationalResults(self.recipe, directory=self.directory)

        if self.recipe['flavor'] == 'gaussian':
            look_for = ['*.fchk']
        else:
            look_for = ['*.tar.gz']
            if any(a[0] == 'G' for a in self.recipe.bases()):
                look_for.append('*.out')

        for directory in directories:
            if not os.path.isdir(directory):
                raise BadCooking('{} is no directory!'.format(directory))

            for l in look_for:
                for i in glob.glob('{}/{}'.format(directory, l)):
                    if verbosity_level >= 1:
                        out.write('* cooking with {} ... '.format(i))

                    with open(i) as f:
                        try:
                            fx = helpers.open_chemistry_file(f)
                            self.cook_from_file(fx, i, storage)
                            if verbosity_level >= 1:
                                out.write('ok\n')
                        except helpers.ProbablyNotAChemistryFile:
                            if verbosity_level >= 1:
                                out.write('skipped\n')
                            continue

        return storage

    def cook_from_file(self, f, name, storage):
        """

        :param f: file
        :type f: qcip_tools.chemistry_files.ChemistryFile
        :param name: path to the file
        :type name: str
        :param storage: storage object
        :type storage: nachos.core.files.ComputationalResults
        :return:
        """

        def almost_the_same(a, b, threshold=1e-3):
            return math.fabs(a - b) < threshold

        # catch field
        if self.recipe['type'] == 'F':
            try:
                real_fields = f.property('n:input_electric_field')[1:4]
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
            return None  # not part of this calculation

        level = next(a[1] for a in self.fields_needed_by_recipe if a[0] == fields)
        derivatives_per_level = self.recipe.bases(level_min=level)
        derivatives_in_level = [a[0] for a in derivatives_per_level]

        if derivatives.Derivative() in derivatives_in_level and f.file_type != 'DALTON_LOG':
            try:
                energies = f.property('computed_energies')
                storage.add_result(
                    fields,
                    '',
                    derivatives.Tensor('', components=numpy.array((energies['total'],))),
                    allow_replace=True)
            except (PropertyNotPresent, PropertyNotDefined):
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

    @staticmethod
    def real_fields_to_fields(real_field, min_field, ratio):
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
    def real_fields_from_geometry(geometry, deformed_geometry, threshold=1e-4):
        if len(geometry) != len(deformed_geometry) or str(geometry) != str(deformed_geometry):
            raise ValueError('geometries does not contain the same number of atoms')

        real_fields = [0] * (3 * len(geometry))

        for index, a in enumerate(deformed_geometry):
            for i in range(3):
                diff_value = (a.position[i] - geometry[index].position[i]) / quantities.AuToAngstrom
                real_fields[index * 3 + i] = 0 if math.fabs(diff_value) < threshold else diff_value

        return real_fields
