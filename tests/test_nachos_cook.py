import os
import random
import subprocess

from qcip_tools.chemistry_files import gaussian, dalton, xyz

from tests import NachosTestCase
from nachos.core import files, cooking, preparing


class CookTestCase(NachosTestCase):

    def setUp(self):
        self.zip_F = self.copy_to_temporary_directory('numdiff_F.zip')
        self.zip_G = self.copy_to_temporary_directory('numdiff_G.zip')
        self.zip_G_dalton = self.copy_to_temporary_directory('numdiff_G_dalton.zip')
        self.working_directory = self.setup_temporary_directory()

    def tearDown(self):
        super().tearDown()
        pass

    def test_fields_from_deformed_geometry(self):
        """Test that the code is able to get back the fields from a deformed geometry"""

        min_field = 0.01
        ratio = 2.

        path = self.copy_to_temporary_directory('water.xyz')
        fx = xyz.File()

        with open(path) as f:
            fx.read(f)

        geometry = fx.molecule
        dof = 3 * len(geometry)

        for _ in range(5):
            fields = [0] * dof
            n = random.randrange(1, 3)
            for i in range(n):
                fields[random.randrange(0, dof)] = random.randrange(-5, 5)
            real_fields = preparing.Preparer.real_fields(fields, min_field, ratio)
            deformed_geometry = preparing.Preparer.deform_geometry(geometry, real_fields)

            self.assertArrayAlmostEqual(
                real_fields, cooking.Cooker.real_fields_from_geometry(geometry, deformed_geometry))

            self.assertArrayAlmostEqual(
                fields, cooking.Cooker.real_fields_to_fields(real_fields, min_field, ratio))

    def test_cook_F_gaussian(self):
        self.unzip_it(self.zip_F, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_F')
        path = os.path.join(directory, 'nachos_recipe.yml')

        r = files.Recipe(directory=directory)

        with open(path) as f:
            r.read(f)

        fields = preparing.fields_needed_by_recipe(r)

        c = cooking.Cooker(r, directory)
        storage = c.cook()

        # write and read
        self.assertEqual(storage.check(), ([], []))
        storage.write('nachos_data.h5')
        s = files.ComputationalResults(r, directory=directory)
        s.read('nachos_data.h5')
        self.assertEqual(s.check(), ([], []))

        # check data
        for _ in range(10):
            n = random.randrange(1, len(fields) + 1)
            fields_n, level = fields[n - 1]
            t_fields = tuple(fields_n)
            path = os.path.join(directory, r['name'] + '_{:04d}.fchk').format(n)
            self.assertTrue(os.path.exists(path), msg=path)
            with open(path) as f:
                fx = gaussian.FCHK()
                fx.read(f)
                results = s.results[t_fields]

                self.assertAlmostEqual(fx.property('computed_energies')['total'], results[''][0])

                if level <= 2:
                    electrical_derivatives = fx.property('electrical_derivatives')
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['F']['static'].components, results['F']['static'])

                    if level <= 1:
                        self.assertArrayAlmostEqual(
                            electrical_derivatives['FF']['static'].components, results['FF']['static'])
                        self.assertArrayAlmostEqual(
                            electrical_derivatives['FD'][0.0428227067].components, results['FD']['1064nm'])

    def test_cook_G_gaussian(self):
        self.unzip_it(self.zip_G, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_G')

        path = os.path.join(directory, 'nachos_recipe.yml')

        r = files.Recipe(directory=directory)

        with open(path) as f:
            r.read(f)

        fields = preparing.fields_needed_by_recipe(r)

        c = cooking.Cooker(r, directory)
        storage = c.cook()

        # write and read
        self.assertEqual(storage.check(), ([], []))
        storage.write('nachos_data.h5')
        s = files.ComputationalResults(r, directory=directory)
        s.read('nachos_data.h5')
        self.assertEqual(s.check(), ([], []))

        for _ in range(10):
            n = random.randrange(1, len(fields) + 1)
            fields_n, level = fields[n - 1]

            path = os.path.join(directory, r['name'] + '_{:04d}{}.fchk').format(n, 'a' if level < 2 else '')
            self.assertTrue(os.path.exists(path), msg=path)
            with open(path) as f:
                fx = gaussian.FCHK()
                fx.read(f)

                t_fields = tuple(cooking.Cooker.real_fields_to_fields(
                    cooking.Cooker.real_fields_from_geometry(r.geometry, fx.molecule), r['min_field'], r['ratio']))

                results = s.results[t_fields]
                self.assertAlmostEqual(fx.property('computed_energies')['total'], results[''][0])

                if level <= 1:
                    electrical_derivatives = fx.property('electrical_derivatives')
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['F']['static'].components, results['F']['static'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FF']['static'].components, results['FF']['static'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FD'][0.0428227067].components, results['FD']['1064nm'])

            path = os.path.join(directory, r['name'] + '_{:04d}{}.fchk').format(n, 'b' if level < 2 else '')
            self.assertTrue(os.path.exists(path), msg=path)
            with open(path) as f:
                fx = gaussian.FCHK()
                fx.read(f)

                if level <= 1:
                    geometrical_derivatives = fx.property('geometrical_derivatives')
                    self.assertArrayAlmostEqual(
                        geometrical_derivatives['G'].components,
                        results['G'])

    def test_cook_G_dalton(self):
        self.unzip_it(self.zip_G_dalton, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_G_dalton')

        path = os.path.join(directory, 'nachos_recipe.yml')

        r = files.Recipe(directory=directory)

        with open(path) as f:
            r.read(f)

        fields = preparing.fields_needed_by_recipe(r)

        c = cooking.Cooker(r, directory)
        storage = c.cook()

        # write and read
        self.assertEqual(storage.check(), ([], []))
        storage.write('nachos_data.h5')
        s = files.ComputationalResults(r, directory=directory)
        s.read('nachos_data.h5')
        self.assertEqual(s.check(), ([], []))

        for _ in range(10):
            n = random.randrange(1, len(fields) + 1)
            fields_n, level = fields[n - 1]

            path = os.path.join(directory, 'ND_{:04d}_{}_{:04d}.tar.gz').format(level, r['name'], n)
            self.assertTrue(os.path.exists(path), msg=path)
            with open(path, 'rb') as f:
                fx = dalton.ArchiveOutput()
                fx.read(f)

                t_fields = tuple(cooking.Cooker.real_fields_to_fields(
                    cooking.Cooker.real_fields_from_geometry(r.geometry, fx.molecule), r['min_field'], r['ratio']))

                results = s.results[t_fields]
                self.assertAlmostEqual(fx.property('computed_energies')['total'], results[''][0])

                if level <= 1:
                    electrical_derivatives = fx.property('electrical_derivatives')
                    fr = 0.04282271
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['F']['static'].components, results['F']['static'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FF']['static'].components, results['FF']['static'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FD'][fr].components, results['FD']['1064nm'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FFF']['static'].components, results['FFF']['static'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FDF'][fr].components, results['FDF']['1064nm'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FDD'][fr].components, results['FDD']['1064nm'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FFFF']['static'].components, results['FFFF']['static'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FDFF'][fr].components, results['FDFF']['1064nm'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FDDF'][fr].components, results['FDDF']['1064nm'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FDDd'][fr].components, results['FDDd']['1064nm'])
                    self.assertArrayAlmostEqual(
                        electrical_derivatives['FDDD'][fr].components, results['FDDD']['1064nm'])

    def test_nachos_cook(self):
        """Test the cooker program"""

        self.unzip_it(self.zip_F, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_F')
        recipe_path = os.path.join(directory, 'nachos_recipe.yml')

        h5_path = os.path.join(directory, 'xxx.h5')
        self.assertFalse(os.path.exists(h5_path))

        # process without copy
        process = self.run_python_script(
            'nachos/nachos_cook.py', ['-r', recipe_path, '-o', h5_path],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertEqual(len(stdout_t), 0, msg=stdout_t.decode())

        self.assertTrue(os.path.exists(h5_path))

        # check h5 file
        r = files.Recipe(directory=os.path.dirname(recipe_path))
        with open(recipe_path) as f:
            r.read(f)

        s = files.ComputationalResults(r, directory=os.path.dirname(recipe_path))
        s.read(h5_path)
        self.assertEqual(s.check(), ([], []))
