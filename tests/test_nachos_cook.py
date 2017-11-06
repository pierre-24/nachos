import os

from tests import NachosTestCase
from nachos.core import files, cooking


class CookTestCase(NachosTestCase):

    def setUp(self):
        self.zip_F = self.copy_to_temporary_directory('numdiff_F.zip')
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

        cooking.Cooker(r, directory)
