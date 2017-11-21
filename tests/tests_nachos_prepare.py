import os
import random
import subprocess
import glob

from qcip_tools import quantities, derivatives_e, numerical_differentiation
from qcip_tools.chemistry_files import gaussian, dalton

from tests import NachosTestCase
from nachos.core import files, preparing


class PrepareTestCase(NachosTestCase):

    def setUp(self):
        self.geometry = self.copy_to_temporary_directory('water.xyz')
        self.basis_set = self.copy_to_temporary_directory('sto-3g.gbs')
        self.custom_recipe = self.copy_to_temporary_directory('nachos_recipe.yml')
        self.working_directory = self.setup_temporary_directory()

    def tearDown(self):
        super().tearDown()
        pass

    def test_fields_needed(self):
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
        self.assertEqual(len(preparing.fields_needed_by_recipe(r)), 137)
        self.assertEqual(preparing.fields_needed_by_recipe(r)[0], ([0, 0, 0], 1))  # zero field is the first one

        r['k_max'] = 3
        self.assertEqual(len(preparing.fields_needed_by_recipe(r)), 85)

        r['accuracy_level'] = 0
        self.assertEqual(len(preparing.fields_needed_by_recipe(r)), 59)

        r['differentiation'] = {1: ['energy']}
        r._update(r.recipe)
        self.assertEqual(len(preparing.fields_needed_by_recipe(r)), 19)

        r['type'] = 'G'
        r['differentiation'] = {2: ['energy']}
        self.assertEqual(len(preparing.fields_needed_by_recipe(r)), 3 * len(r.geometry) * 2 * r['k_max'] + 1)
        self.assertEqual(preparing.fields_needed_by_recipe(r)[0], ([0] * 3 * len(r.geometry), 1))  # zero field

        r['k_max'] = 5
        r._update(r.recipe)
        self.assertEqual(len(preparing.fields_needed_by_recipe(r)), (3 * len(r.geometry)) ** 2 * 2 * r['k_max'] + 1)

    def test_deform_geometry(self):
        """Test geometry deformation"""

        opt_dict = dict(
            flavor='gaussian',
            type='G',
            method='HF',
            basis_set='gen',
            geometry=self.geometry,
            frequencies=['1064nm'],
            min_field=0.01  # atomic units
        )

        r = files.Recipe(**opt_dict)
        r.check_data()

        zero_fields = [0] * r.dof

        # nothing is deformed
        deformed = preparing.Preparer.deform_geometry(r.geometry, zero_fields)
        for i, a in enumerate(deformed):
            self.assertArraysAlmostEqual(a.position, r.geometry[i].position)

        # simple deformation:
        fields = zero_fields.copy()
        fields[0] = 1  # x of first atom
        deformed = preparing.Preparer.deform_geometry(
            r.geometry, numerical_differentiation.real_fields(fields, r['min_field'], r['ratio']))

        self.assertNotEqual(r.geometry[0].position[0], deformed[0].position[0])
        self.assertEqual(r.geometry[0].position[1], deformed[0].position[1])  # only the first coordinate was modified
        self.assertEqual(r.geometry[0].position[2], deformed[0].position[2])

        for i, a in enumerate(deformed[1:]):
            self.assertArraysAlmostEqual(a.position, r.geometry[i + 1].position)  # nothing else was deformed !

        self.assertEqual(r.geometry[0].position[0] + r['min_field'] * quantities.AuToAngstrom, deformed[0].position[0])

        # double deformation:
        fields = zero_fields.copy()
        fields[0] = 3
        fields[4] = -2  # y of second atom
        deformed = preparing.Preparer.deform_geometry(
            r.geometry, numerical_differentiation.real_fields(fields, r['min_field'], r['ratio']))

        self.assertNotEqual(r.geometry[0].position[0], deformed[0].position[0])
        self.assertNotEqual(r.geometry[1].position[1], deformed[1].position[1])

        self.assertEqual(
            r.geometry[0].position[0] + r['min_field'] * r['ratio'] ** 2 * quantities.AuToAngstrom,
            deformed[0].position[0])

        self.assertEqual(
            r.geometry[1].position[1] - r['min_field'] * r['ratio'] * quantities.AuToAngstrom,
            deformed[1].position[1])

    def test_preparer_for_F(self):
        """Test the preparer class"""

        name = 'water_test'
        min_field = .0004

        differentiation = {
            3: ['energy'],
            2: ['F'],
            1: ['FF', 'FD']
        }

        opt_dict = dict(
            flavor='gaussian',
            type='F',
            method='HF',
            basis_set='gen',
            geometry=self.geometry,
            differentiation=differentiation,
            frequencies=['1064nm'],
            name=name,
            min_field=min_field
        )

        flavor_extra = {
            'memory': '3Gb',
            'procs': 4,
            'gen_basis': self.basis_set
        }

        r = files.Recipe(flavor_extra=flavor_extra, **opt_dict)
        r.check_data()

        fields = preparing.fields_needed_by_recipe(r)

        preparer = preparing.Preparer(recipe=r, directory=self.working_directory)
        preparer.prepare()

        # test for base
        path = os.path.join(self.working_directory, name + '_0001.com')
        self.assertTrue(os.path.exists(path))
        with open(path) as f:
            fi = gaussian.Input()
            fi.read(f)

            self.assertEqual(len(fi.other_blocks), 3)
            self.assertAlmostEqual(float(fi.other_blocks[0][0]), derivatives_e.convert_frequency_from_string('1064nm'))
            self.assertEqual(fi.other_blocks[1][0], 'O     0')
            self.assertEqual([float(a) for a in fi.other_blocks[2][0].split()], [.0, .0, .0])  # zero field

        for _ in range(5):  # 5 random tests that files contains what they should!
            n = random.randrange(1, len(fields) + 1)
            path = os.path.join(self.working_directory, name + '_{:04d}.com').format(n)
            self.assertTrue(os.path.exists(path), msg=path)
            with open(path) as f:
                fi = gaussian.Input()
                fi.read(f)

                fields_n, level = fields[n - 1]

                self.assertEqual(len(fi.other_blocks), 3 if level < 2 else 2)
                if level < 2:
                    self.assertAlmostEqual(
                        float(fi.other_blocks[0][0]), derivatives_e.convert_frequency_from_string('1064nm'))

                self.assertEqual(fi.other_blocks[-2][0], 'O     0')
                self.assertEqual(
                    [float(a) for a in fi.other_blocks[-1][0].split()],
                    numerical_differentiation.real_fields(fields_n, min_field, 2.))

    def test_preparer_for_G(self):
        """Test the preparer class"""

        name = 'water_test'
        min_field = .01

        differentiation = {
            2: ['energy'],
            1: ['G', 'F', 'FF', 'FD']
        }

        opt_dict = dict(
            flavor='gaussian',
            type='G',
            method='HF',
            basis_set='gen',
            geometry=self.geometry,
            differentiation=differentiation,
            frequencies=['1064nm'],
            name=name,
            min_field=min_field,
            k_max=3
        )

        flavor_extra = {
            'memory': '3Gb',
            'procs': 4,
            'gen_basis': self.basis_set
        }

        r = files.Recipe(flavor_extra=flavor_extra, **opt_dict)
        r.check_data()

        fields = preparing.fields_needed_by_recipe(r)

        recipe_path = os.path.join(self.working_directory, 'nachos_recipe.yml')
        with open(recipe_path, 'w') as f:
            r.write(f)

        preparer = preparing.Preparer(recipe=r, directory=self.working_directory)
        preparer.prepare()

        path = os.path.join(self.working_directory, name + '_0001a.com')
        self.assertTrue(os.path.exists(path))

        # test base files
        with open(path) as f:
            fi = gaussian.Input()
            fi.read(f)

            self.assertTrue(any(['polar' in a for a in fi.input_card]))
            self.assertEqual(len(fi.other_blocks), 2)
            self.assertAlmostEqual(float(fi.other_blocks[0][0]), derivatives_e.convert_frequency_from_string('1064nm'))
            self.assertEqual(fi.other_blocks[1][0], 'O     0')

            for i, a in enumerate(fi.molecule):
                self.assertArraysAlmostEqual(r.geometry[i].position, a.position)  # geometry not modified, it is base

        path = os.path.join(self.working_directory, name + '_0001b.com')
        self.assertTrue(os.path.exists(path))
        with open(path) as f:
            fi = gaussian.Input()
            fi.read(f)

            self.assertFalse(any(['polar' in a for a in fi.input_card]))
            self.assertEqual(len(fi.other_blocks), 1)
            self.assertEqual(fi.other_blocks[0][0], 'O     0')

            for i, a in enumerate(fi.molecule):
                self.assertArraysAlmostEqual(r.geometry[i].position, a.position)  # geometry not modified, it is base

        for _ in range(5):  # 5 random tests that files contains what they should!
            n = random.randrange(1, len(fields) + 1)
            fields_n, level = fields[n - 1]
            path = os.path.join(self.working_directory, name + '_{:04d}{}.com').format(n, 'a' if level < 2 else '')
            self.assertTrue(os.path.exists(path), msg=path)
            with open(path) as f:
                fi = gaussian.Input()
                fi.read(f)

                if level < 2:
                    self.assertTrue(any(['polar' in a for a in fi.input_card]))
                    self.assertAlmostEqual(
                        float(fi.other_blocks[0][0]), derivatives_e.convert_frequency_from_string('1064nm'))
                else:
                    self.assertFalse(any(['polar' in a for a in fi.input_card]))

                self.assertEqual(fi.other_blocks[-1][0], 'O     0')

                deformed = preparing.Preparer.deform_geometry(
                    r.geometry, numerical_differentiation.real_fields(fields_n, min_field, 2))

                for i, a in enumerate(fi.molecule):
                    self.assertArraysAlmostEqual(deformed[i].position, a.position)

    def test_preparer_for_G_dalton(self):
        """Test the preparer class, for dalton"""

        name = 'water_test'
        min_field = .01

        differentiation = {
            2: ['energy'],
            1: ['G', 'F', 'FF', 'FD', 'FFF', 'FDF', 'FDD', 'FFFF', 'FDFF', 'FDDF', 'FDDd', 'FDDD']
        }

        opt_dict = dict(
            flavor='dalton',
            type='G',
            method='CCS',
            basis_set='STO-3G',
            geometry=self.geometry,
            differentiation=differentiation,
            frequencies=['1064nm'],
            name=name,
            min_field=min_field,
            k_max=3
        )

        r = files.Recipe(**opt_dict)
        r.check_data()

        fields = preparing.fields_needed_by_recipe(r)

        recipe_path = os.path.join(self.working_directory, 'nachos_recipe.yml')
        with open(recipe_path, 'w') as f:
            r.write(f)

        preparer = preparing.Preparer(recipe=r, directory=self.working_directory)
        preparer.prepare()

        path = os.path.join(self.working_directory, name + '_0001.mol')
        self.assertTrue(os.path.exists(path))

        # test base files
        with open(path) as f:
            fi = dalton.MoleculeInput()
            fi.read(f)

            self.assertEqual(fi.basis_set, 'STO-3G')

            for i, a in enumerate(fi.molecule):
                self.assertArraysAlmostEqual(r.geometry[i].position, a.position)  # geometry not modified, it is base

        for _ in range(5):  # 5 random tests that files contains what they should!
            n = random.randrange(1, len(fields) + 1)
            fields_n, level = fields[n - 1]
            path = os.path.join(self.working_directory, name + '_{:04d}.mol').format(n)
            self.assertTrue(os.path.exists(path), msg=path)
            with open(path) as f:
                fi = dalton.MoleculeInput()
                fi.read(f)

                self.assertEqual(fi.basis_set, 'STO-3G')

                deformed = preparing.Preparer.deform_geometry(
                    r.geometry, numerical_differentiation.real_fields(fields_n, min_field, 2))

                for i, a in enumerate(fi.molecule):
                    self.assertArraysAlmostEqual(deformed[i].position, a.position)

    def test_nachos_prepare(self):
        """Test the preparer program"""

        self.assertEqual(len([a for a in glob.glob(self.working_directory + '/*.com')]), 0)
        self.assertEqual(len([a for a in glob.glob(self.working_directory + '/*.yml')]), 0)
        self.assertEqual(len([a for a in glob.glob(self.working_directory + '/*.xyz')]), 0)
        self.assertEqual(len([a for a in glob.glob(self.working_directory + '/*.gbs')]), 0)

        # process without copy
        process = self.run_python_script(
            'nachos/prepare.py',
            ['-r', self.custom_recipe, '-d', self.working_directory],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertNotEqual(len(stdout_t), 0)

        self.assertEqual(len([a for a in glob.glob(self.working_directory + '/*.com')]), 542)
        self.assertEqual(len([a for a in glob.glob(self.working_directory + '/*.yml')]), 0)
        self.assertEqual(len([a for a in glob.glob(self.working_directory + '/*.xyz')]), 0)
        self.assertEqual(len([a for a in glob.glob(self.working_directory + '/*.gbs')]), 0)

        # process with copy
        other_dir = self.setup_temporary_directory()
        process = self.run_python_script(
            'nachos/prepare.py',
            ['-r', self.custom_recipe, '-d', other_dir, '-c'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertNotEqual(len(stdout_t), 0)

        self.assertEqual(len([a for a in glob.glob(other_dir + '/*.com')]), 542)
        self.assertEqual(len([a for a in glob.glob(other_dir + '/*.yml')]), 1)
        self.assertEqual(len([a for a in glob.glob(other_dir + '/*.xyz')]), 1)
        self.assertEqual(len([a for a in glob.glob(other_dir + '/*.gbs')]), 1)
