import os

from tests import NachosTestCase
from nachos.core import files


class FilesTestCase(NachosTestCase):

    def setUp(self):
        self.geometry = self.copy_to_temporary_directory('water.xyz')

    def test_recipe(self):
        """Test the behavior of a recipe"""

        opt_dict = dict(
            flavor='gaussian',
            type='electric',
            method='HF',
            basis_set='STO-3G',
            geometry=self.geometry)

        flavor_extra = {
            'memory': '3Gb',
            'procs': 4
        }

        r = files.Recipe(flavor_extra=flavor_extra, **opt_dict)

        for key in opt_dict:
            self.assertEqual(r.recipe[key], opt_dict[key])

        for key in flavor_extra:
            self.assertEqual(r.recipe['flavor_extra'][key], flavor_extra[key])

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

        with self.assertRaises(files.BadRecipe):
            opt_dict_p = opt_dict.copy()
            opt_dict_p['bases'] = ['energy', 'F', 'FD']
            rx = files.Recipe(**opt_dict_p)
            rx.check_data()

        with self.assertRaises(files.BadRecipe):
            rx = files.Recipe(flavor_extra={'x': 'y'}, **opt_dict)
            rx.check_data()
