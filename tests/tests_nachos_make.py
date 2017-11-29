import subprocess
import os

from nachos.core import files

from tests import NachosTestCase


class MakeTestCase(NachosTestCase):

    def setUp(self):
        self.geometry = self.copy_to_temporary_directory('water.xyz')
        self.basis_set = self.copy_to_temporary_directory('sto-3g.gbs')
        self.working_directory = self.setup_temporary_directory()

    def tearDown(self):
        super().tearDown()
        pass

    def test_nachos_make(self):
        """Test the maker program"""

        opt_dict = dict(
            flavor='gaussian',
            type='F',
            method='HF',
            basis_set='STO-3G',
            geometry=self.geometry,
            differentiation={3: ['energy', 'FD']},
            frequencies=['1064nm'],
            accuracy_level=1,
            k_max=5,
            ratio=2,
            min_field=.001
        )

        r = files.Recipe(**opt_dict)

        recipe_path = os.path.join(self.working_directory, 'recipe.yml')

        opt_list = [
            '-N',
            '--flavor', r['flavor'],
            '--type', r['type'],
            '--method', r['method'],
            '--basis-set', r['basis_set'],
            '--geometry', self.geometry,
            '--differentiation', 'energy:3;FD:3',
            '--frequencies', '1064nm',
            '--name', r['name'],
            '--k-max', str(r['k_max']),
            '--ratio', str(r['ratio']),
            '--min-field', str(r['min_field']),
            '--flavor-extra', '',
            '-o', recipe_path
        ]

        process = self.run_python_script(
            'nachos/make.py',
            opt_list,
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertEqual(len(stdout_t), 0, msg=stdout_t.decode())

        r_generated = files.Recipe(directory=self.working_directory)
        self.assertTrue(os.path.exists(recipe_path))

        with open(recipe_path) as f:
            r_generated.read(f)

        to_check = [
            'flavor', 'type', 'method', 'basis_set', 'geometry', 'frequencies', 'name', 'k_max', 'ratio', 'min_field',
            'differentiation']

        for c in to_check:
            self.assertEqual(r.recipe[c], r_generated.recipe[c], msg=c)
