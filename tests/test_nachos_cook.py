import os
import random

from qcip_tools.chemistry_files import gaussian, dalton

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

    def test_cook_F_gaussian(self):
        self.unzip_it(self.zip_F, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_F')

        path = os.path.join(directory, 'nachos_recipe.yml')

        r = files.Recipe(directory=directory)

        with open(path) as f:
            r.read(f)

        fields = preparing.fields_needed_by_recipe(r)

        c = cooking.Cooker(r, directory)
        c.cook()

        # write and read
        self.assertEqual(c.storage.check(), ([], []))
        c.storage.write('nachos_data.h5')
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
        c.cook()

        # write and read
        self.assertEqual(c.storage.check(), ([], []))
        c.storage.write('nachos_data.h5')
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
        c.cook()

        # write and read
        self.assertEqual(c.storage.check(), ([], []))
        c.storage.write('nachos_data.h5')
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
