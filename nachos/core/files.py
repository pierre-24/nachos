import yaml
import os

from qcip_tools import derivatives
from qcip_tools.chemistry_files import helpers

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

    def __init__(self, **kwargs):

        self.geometry = None

        self.recipe = {}
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

        if not os.path.exists(self['geometry']):
            raise BadRecipe('file "{}" containing the geometry does not exists'.format(self['geometry']))

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
        if 'geometry' in kw and os.path.exists(kw['geometry']):
            try:
                with open(kw['geometry']) as f:
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


class ComputationalResults:
    """A class to store all the results obtained from the baking process

    :param recipe: a recipe
    :type recipe: nachos.core.files.Recipe
    """

    def __init__(self, recipe):
        self.recipe = recipe
