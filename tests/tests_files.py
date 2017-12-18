import os

from tests import NachosTestCase
from nachos.core import files


class FilesTestCase(NachosTestCase):

    def setUp(self):
        self.geometry = self.copy_to_temporary_directory('water.xyz')

    def test_recipe(self):
        """Test the behavior of a recipe"""

        differentiation = {
            3: ['energy'],
            2: ['F'],
            1: ['FF']
        }

        opt_dict = dict(
            flavor='gaussian',
            type='F',
            method='HF',
            basis_set='STO-3G',
            geometry=self.geometry,
            differentiation=differentiation
        )

        flavor_extra = {
            'memory': '3Gb',
            'procs': 4
        }

        r = files.Recipe(flavor_extra=flavor_extra, **opt_dict)

        for key in opt_dict:
            self.assertEqual(r.recipe[key], opt_dict[key])

        for key in flavor_extra:
            self.assertEqual(r.recipe['flavor_extra'][key], flavor_extra[key])

        bases = r.bases()
        for basis, level in bases:
            representation = basis.representation()
            if representation == '':
                representation = 'energy'
            self.assertIn(representation, differentiation[level])

        # test write
        recipe_file = os.path.join(self.temporary_directory, 'nachos_recipe.yaml')
        with open(recipe_file, 'w') as f:
            r.write(f)

        rp = files.Recipe()
        self.assertEqual(rp['flavor_extra'], {})

        with open(recipe_file) as f:
            rp.read(f)
            rp.check_data()

        for key in opt_dict:
            self.assertEqual(rp.recipe[key], opt_dict[key])

        for key in flavor_extra:
            self.assertEqual(rp.recipe['flavor_extra'][key], flavor_extra[key])

        # some exceptions
        with self.assertRaises(files.BadRecipe):
            rx = files.Recipe()
            rx.check_data()

        with self.assertRaises(files.BadRecipe):
            rx = files.Recipe(flavor='x')
            rx.check_data()

        with self.assertRaises(files.BadRecipe):
            rx = files.Recipe(accuracy_level=10000, **opt_dict)
            rx.check_data()

        with self.assertRaises(files.BadRecipe):  # differentiation level should be an int
            opt_dict_p = opt_dict.copy()
            opt_dict_p['differentiation'] = {'1': ['energy', 'F']}
            rx = files.Recipe(**opt_dict_p)
            rx.check_data()

        with self.assertRaises(files.BadRecipe):  # differentiation level should be larger than 1
            opt_dict_p = opt_dict.copy()
            opt_dict_p['differentiation'] = {0: ['energy', 'F']}
            rx = files.Recipe(**opt_dict_p)
            rx.check_data()

        with self.assertRaises(files.BadRecipe):  # no frequency
            opt_dict_p = opt_dict.copy()
            opt_dict_p['differentiation'] = {1: ['energy', 'F', 'FD']}
            rx = files.Recipe(**opt_dict_p)
            rx.check_data()

        with self.assertRaises(files.BadRecipe):  # unknown option for flavor
            rx = files.Recipe(flavor_extra={'x': 'y'}, **opt_dict)
            rx.check_data()

        with self.assertRaises(files.BadRecipe):  # polarizability is not available for MP3
            opt_dict_p = opt_dict.copy()
            opt_dict_p['method'] = 'MP3'
            opt_dict_p['differentiation'] = {1: ['energy', 'F', 'FF']}
            rx = files.Recipe(**opt_dict_p)
            rx.check_data()
