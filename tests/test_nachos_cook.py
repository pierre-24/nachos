import os

from tests import NachosTestCase
from nachos.core import files, cooking


class CookTestCase(NachosTestCase):

    def setUp(self):
        self.zip_F = self.copy_to_temporary_directory('numdiff_F.zip')
        self.zip_G = self.copy_to_temporary_directory('numdiff_G.zip')
        self.zip_G_dalton = self.copy_to_temporary_directory('numdiff_G_dalton.zip')
        self.working_directory = self.setup_temporary_directory()

    def tearDown(self):
        super().tearDown()
        pass

    def test_cook_F_gaussian(self):
        self.unzip_it(self.zip_F, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_F')

        path = os.path.join(directory, 'nachos_recipe.yml')

        r = files.Recipe(directory=directory)

        with open(path) as f:
            r.read(f)

        c = cooking.Cooker(r, directory)
        c.cook()

        # write and read
        self.assertEqual(c.storage.check(), ([], []))
        c.storage.write('nachos_data.h5')
        s = files.ComputationalResults(r, directory=directory)
        s.read('nachos_data.h5')
        self.assertEqual(s.check(), ([], []))

    def test_cook_G_gaussian(self):
        self.unzip_it(self.zip_G, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_G')

        path = os.path.join(directory, 'nachos_recipe.yml')

        r = files.Recipe(directory=directory)

        with open(path) as f:
            r.read(f)

        c = cooking.Cooker(r, directory)
        c.cook()

        # write and read
        self.assertEqual(c.storage.check(), ([], []))
        c.storage.write('nachos_data.h5')
        s = files.ComputationalResults(r, directory=directory)
        s.read('nachos_data.h5')
        self.assertEqual(s.check(), ([], []))

    def test_cook_G_dalton(self):
        self.unzip_it(self.zip_G_dalton, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_G_dalton')

        path = os.path.join(directory, 'nachos_recipe.yml')

        r = files.Recipe(directory=directory)

        with open(path) as f:
            r.read(f)

        c = cooking.Cooker(r, directory)
        c.cook()

        # write and read
        self.assertEqual(c.storage.check(), ([], []))
        c.storage.write('nachos_data.h5')
        s = files.ComputationalResults(r, directory=directory)
        s.read('nachos_data.h5')
        self.assertEqual(s.check(), ([], []))
