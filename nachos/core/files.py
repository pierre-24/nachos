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
    'bases',
    'max_differentiation',
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
    'bases': ['energy'],
    'max_differentiation': 3,
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
        self.recipe.update(**DEFAULT_RECIPE)
        self.recipe.update(**kwargs)
        self.__update(kwargs)

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
        if self['method'] not in CONFIG[self['flavor']]['methods']:
            raise BadRecipe('Calculation with method {} not available'.format(self['method']))

        for basis in self['bases']:
            if basis not in config_for_flavor['bases']:
                raise BadRecipe('Differentiation of {} not available'.format(basis))

            if 'D' in basis and 'frequencies' not in self.recipe:
                raise BadRecipe('Field dependant basis requested ({}) but no "frequencies" field'.format(basis))

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
        self.__update(up)
        self.check_data()

    def __update(self, kw):
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
            except helpers.ProbablyNotAChemistryFile:
                pass

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

    def bases(self):
        """Get basis in the form of Derivatives objects

        :rtype: list
        """

        dof = 3 * len(self.geometry)
        return [derivatives.Derivative(basis if basis != 'energy' else '', spacial_dof=dof) for basis in self['bases']]

    def derivatives(self):
        """Get the different derivatives performed by the recipe

        :rtype: list
        """

        diff_repr = 'F' if self['type'] == 'electric' else 'G'
        d = []

        for i in range(1, self['max_differentiation'] + 1):
            d.append(derivatives.Derivative(diff_repr * i, spacial_dof=3 * len(self.geometry)))

        return d


class ComputationalResults:
    """A class to store all the results obtained from the baking process

    :param recipe: a recipe
    :type recipe: nachos.core.files.Recipe
    """

    def __init__(self, recipe):
        self.recipe = recipe
