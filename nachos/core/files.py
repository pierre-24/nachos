import yaml
import os
import h5py
import numpy

from nachos.core import preparing

from qcip_tools import derivatives
from qcip_tools.chemistry_files import helpers, chemistry_datafile

from nachos.core import CONFIG

RECIPE_MANDATORY_FIELDS = [
    'flavor',
    'type',
    'method',
    'basis_set',
    'geometry',

    # fields for which there is default:
    'name',
    'min_field',
    'ratio',
    'k_max',
    'differentiation',
    'accuracy_level',  # -1 = T-REX compatibility, 0 = default, 1 = more accurate for order 3
    'flavor_extra'
]

RECIPE_OPTIONAL_FIELDS = [
    'frequencies',
]

DEFAULT_RECIPE = {
    'name': 'ND',
    'min_field': 0.0004,
    'ratio': 2,
    'k_max': 5,
    'differentiation': {2: ['energy']},
    'accuracy_level': 1,
    'flavor_extra': {},
}


class BadRecipe(Exception):
    pass


class Recipe:
    """Class that handle the parameters to perform numerical differentiation"""

    def __init__(self, directory='.', **kwargs):

        self.geometry = None

        if not os.path.isdir(directory):
            raise BadRecipe('{} is not a directory'.format(directory))

        self.recipe = {}
        self.directory = directory
        self.max_differentiation = 1
        self.dof = 0
        self.recipe.update(**DEFAULT_RECIPE)
        self.recipe.update(**kwargs)
        self._update(kwargs)

    def check_data(self):
        """check if data are coherent

        :raise BadRecipe: if there is something wrong ;)
        """

        for key in self.recipe:
            if key not in RECIPE_MANDATORY_FIELDS and key not in RECIPE_OPTIONAL_FIELDS:
                raise BadRecipe('field "{}" is not allowed'.format(key))

        for mandatory_field in RECIPE_MANDATORY_FIELDS:
            if mandatory_field not in self.recipe:
                raise BadRecipe('field {} is missing'.format(mandatory_field))

        if self['accuracy_level'] not in [-1, 0, 1]:
            raise BadRecipe('Accuracy level {} does not exists'.format(self['accuracy_level']))

        path = os.path.join(self.directory, self['geometry'])

        if not os.path.exists(path):
            raise BadRecipe('file "{}" containing the geometry does not exists'.format(path))

        if self.geometry is None:
            raise BadRecipe('geometry is not readable')

        if self['flavor'] not in CONFIG:
            raise BadRecipe('Flavor {} not allowed'.format(self['flavor']))

        config_for_flavor = CONFIG[self['flavor']]

        if self['type'] not in CONFIG[self['flavor']]['types']:
            raise BadRecipe('{}al derivatives not available for flavor "{}"'.format(self['type'], self['flavor']))
        if self['method'] not in [a[0] for a in CONFIG[self['flavor']]['methods']]:
            raise BadRecipe('Calculation with method {} not available'.format(self['method']))

        max_differentiation_for_method = next(
            a for a in CONFIG[self['flavor']]['methods'] if a[0] == self['method'])[1]

        for level in self['differentiation']:
            if type(level) is not int:
                raise BadRecipe('level of differentiation should be an int')
            if level < 1:
                raise BadRecipe('level of differentiation should be > 1')
            if type(self['differentiation'][level]) is not list:
                raise BadRecipe('bases should be given as a list!')

            for basis in self['differentiation'][level]:
                if basis not in config_for_flavor['bases']:
                    raise BadRecipe('Differentiation of {} not available'.format(basis))

                if 'D' in basis and 'frequencies' not in self.recipe:
                    raise BadRecipe('Field dependant basis requested ({}) but no "frequencies" field'.format(basis))

                num_G = basis.count('G')
                num_F = basis.count('F') + basis.count('D')

                if num_G > max_differentiation_for_method['G']:
                    raise BadRecipe('"{}" for method {} is not possible'.format(basis, self['method']))

                if num_F > max_differentiation_for_method['F']:
                    raise BadRecipe('"{}" for method {} is not possible'.format(basis, self['method']))

        if 'flavor_extra' in self:
            extra = self['flavor_extra']

            for key in extra:
                if key not in config_for_flavor['default_for_extra_fields']:
                    raise BadRecipe(
                        'field {} not allowed in "flavor_extra" for flavor "{}"'.format(key, self['flavor']))

        if self['ratio'] <= 1.0:
            raise BadRecipe('ratio must be > 1.0')

        if self['k_max'] < 1:
            raise BadRecipe('k_max should be > 1')

    def read(self, fp):
        """Load a given YAML file and fill the Recipe with it.

        :param fp: valid file descriptor (in read mode)
        :type fp: file
        :raise BadRecipe: if parameters are not allowed
        """

        up = yaml.load(fp)
        self.recipe.update(up)
        self._update(up)
        self.check_data()

    def _update(self, kw):
        """Set up the default options depending on the flavor, then update their value with what is inside ``kw``."""

        # flavor
        if 'flavor' in self and self['flavor'] in CONFIG:
            self.recipe['flavor_extra'] = CONFIG[self['flavor']]['default_for_extra_fields'].copy()

            if 'flavor_extra' in kw:
                self.recipe['flavor_extra'].update(kw['flavor_extra'])

        # geometry
        if 'geometry' in kw:
            path = os.path.join(self.directory, kw['geometry'])
            if os.path.exists(path):
                try:
                    with open(path) as f:
                        g = helpers.open_chemistry_file(f)
                        if g.has_property('molecule'):
                            self.geometry = g.property('molecule')
                            self.dof = 3 * len(self.geometry)
                except helpers.ProbablyNotAChemistryFile:
                    pass

        # max diff
        if 'differentiation' in kw:
            self.max_differentiation = 0
            for level in kw['differentiation']:
                if type(level) is int:
                    if level > self.max_differentiation:
                        self.max_differentiation = level

    def write(self, fp):
        """Dump the Recipe into a YAML file

        :param fp: valid file descriptor (in write mode)
        :type fp: file
        """

        self.check_data()
        fp.write(yaml.dump(self.recipe))

    def __getitem__(self, item):
        return self.recipe[item]

    def __setitem__(self, key, value):
        self.recipe[key] = value

    def __contains__(self, item):
        return item in self.recipe

    def bases(self, level_min=-1):
        """Get basis in the form of Derivatives objects, plus the level of differentiation

        :param level_min: minimum level required
        :type level_min: int
        :rtype: list of tuple
        """

        bases = []

        for level in self['differentiation']:
                if level >= level_min:
                    for basis in self['differentiation'][level]:
                        bases.append(
                            (derivatives.Derivative(basis if basis != 'energy' else '', spacial_dof=self.dof), level))

        return bases

    def maximum_derivatives(self):
        """Get the different derivatives performed by the recipe

        :rtype: list
        """

        diff_repr = self['type']
        d = []

        for i in range(1, self.max_differentiation + 1):
            d.append(derivatives.Derivative(diff_repr * i, spacial_dof=self.dof))

        return d


class FieldsNotNeeded(Exception):
    pass


class DerivativeAlreadyDefined(Exception):
    def __init__(self, fields, derivative):
        super().__init__('already defined {}::{}'.format(fields, derivative))


class BadResult(Exception):
    pass


class ComputationalResults:
    """A class to store all the results obtained from the cooking process

    :param recipe: a recipe
    :type recipe: nachos.core.files.Recipe
    """

    file_type = 'NACHOS_CR'
    version = 1

    def __init__(self, recipe, directory='.'):

        if not os.path.isdir(directory):
            raise BadRecipe('{} is not a directory'.format(directory))

        self.recipe = recipe
        self.directory = directory

        self.results = {}
        self.fields_needed_by_recipe = preparing.fields_needed_by_recipe(self.recipe)
        self.fields_needed = [a[0] for a in self.fields_needed_by_recipe]

    def add_result(self, fields, derivative, value, allow_replace=False):
        """Add result for a given derivative in given fields

        :param fields: fields
        :type fields: tuple|list
        :param derivative: derivative
        :type derivative: str
        :param value: value of the derivative
        :type value: dict|qcip_tools.derivatives.Tensor
        :param allow_replace: allow the value to be replaced
        :type allow_replace: bool
        """

        if fields not in self.fields_needed:
            raise FieldsNotNeeded(fields)

        t_fields = tuple(fields)

        if t_fields not in self.results:
            self.results[t_fields] = {}

        if derivative not in self.results[t_fields] or allow_replace:
            self.results[t_fields][derivative] = value
        else:
            raise DerivativeAlreadyDefined(fields, derivative)

    def check(self):
        """Check that, according to the recipe, everything is present

        :rtype: tuple
        """

        missing_fields = []
        missing_derivatives = []

        for fields, level in self.fields_needed_by_recipe:
            t_fields = tuple(fields)
            derivatives_needed = self.recipe.bases(level)

            if t_fields not in self.results:
                missing_fields.append(t_fields)

            else:
                for derivative, lvl in derivatives_needed:
                    if derivative.representation() not in self.results[t_fields]:
                        missing_derivatives.append((t_fields, derivative))

        return missing_fields, missing_derivatives

    def write(self, path):
        """Write in h5 file

        :param path: path to the file, relative to directory
        :type path: str
        """

        dof = 3 * len(self.recipe.geometry)

        with h5py.File(os.path.join(self.directory, path), 'w') as f:
            # first protection: version and type
            dset = f.create_dataset('version', (1,), dtype='i', data=self.version)
            dset.attrs['type'] = self.file_type

            # second protection; some of the information of the recipe
            f.create_dataset(
                'recipe', shape=(6,), dtype='f', data=ComputationalResults.get_recipe_check_data(self.recipe))

            # then the results
            fields_group = f.create_group('fields')

            for fields in self.results:
                subgroup = fields_group.create_group(','.join(str(a) for a in fields))
                chemistry_datafile.ChemistryDataFile.write_derivatives_in_group(subgroup, self.results[fields], dof)

    def read(self, path):
        """Read in h5 file

        :param path: path to the file, relative to directory
        :type path: str
        """

        dof = 3 * len(self.recipe.geometry)

        with h5py.File(os.path.join(self.directory, path), 'r') as f:
            if 'version' not in f or f['version'][0] != self.version:
                raise BadResult('version > 1 (={})'.format(f['version'][0]))

            if 'type' not in f['version'].attrs or f['version'].attrs['type'] != self.file_type:
                raise BadResult('type is incorrect')

            if 'recipe' not in f:
                raise BadResult('no information about the recipe in storage')

            if not numpy.any(f['recipe'][:] - ComputationalResults.get_recipe_check_data(self.recipe)):
                raise BadResult('it is probably a storage from a different recipe that you are using!')

            fields_group = f['/fields']

            for i in fields_group:
                fields = [int(a) for a in i.split(',')]
                if fields in self.fields_needed:
                    t_fields = tuple(fields)
                    self.results[t_fields] = chemistry_datafile.ChemistryDataFile.read_derivatives_from_group(
                        fields_group[i], dof)

    def tensor_element_access(self, fields, min_field, basis, inverse, component, frequency, recipe):
        t_fields = tuple(fields)
        if t_fields not in self.results:
            raise BadResult('fields {} is not available'.format(fields))

        results = self.results[t_fields]
        b_repr = basis.representation()
        if b_repr not in results:
            raise BadResult('derivative {} is not available for {}'.format(b_repr, fields))

        results = results[b_repr]

        if derivatives.is_electrical(basis):
            if frequency not in results:
                raise BadResult('frequency {} is not available for {} of {}'.format(frequency, b_repr, fields))

            results = results[frequency]

        if b_repr != '':
            if len(component) != results.representation.order():
                raise BadResult('shape does not match for {}'.format(b_repr))

        return (-1. if inverse else 1.) * results.components[component]

    @staticmethod
    def get_recipe_check_data(recipe):
        return [
            recipe['min_field'],
            recipe['ratio'],
            recipe['k_max'],
            0. if recipe['type'] == 'F' else 1.,
            len(recipe.geometry),
            sum((a + 1) * b.atomic_number for a, b in enumerate(recipe.geometry))
        ]
