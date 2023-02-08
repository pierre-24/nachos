import os
import random
import subprocess

from qcip_tools import numerical_differentiation
from qcip_tools.chemistry_files import gaussian, dalton, xyz

import nachos.qcip_tools_ext.qchem
from tests import NachosTestCase
from nachos.core import files, cooking, preparing


class CookTestCase(NachosTestCase):

    def setUp(self):
        self.zip_F = 'numdiff_F.zip'
        self.zip_F_qchem = 'numdiff_F_qchem.zip'
        self.zip_F_scs_mp2 = 'numdiff_F_SCS-MP2.zip'
        self.zip_F_b2plyp = 'numdiff_F_b2plyp.zip'
        self.zip_G = 'numdiff_G.zip'
        self.zip_G_dalton = 'numdiff_G_dalton.zip'
        self.working_directory = self.setup_temporary_directory()

    def tearDown(self):
        super().tearDown()

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
            real_fields = numerical_differentiation.real_fields(fields, min_field, ratio)
            deformed_geometry = preparing.Preparer.deform_geometry(geometry, real_fields)

            self.assertArraysAlmostEqual(
                real_fields, cooking.Cooker.real_fields_from_geometry(geometry, deformed_geometry))

            self.assertArraysAlmostEqual(
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
        storage = c.cook([directory])

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

                self.assertAlmostEqual(fx.property('computed_energies')['total'], results[''].components[0])

                if level <= 2:
                    electrical_derivatives = fx.property('electrical_derivatives')
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['F']['static'], results['F']['static'])

                    if level <= 1:
                        self.assertTensorsAlmostEqual(
                            electrical_derivatives['FF']['static'], results['FF']['static'])

                        self.assertTensorsAlmostEqual(
                            electrical_derivatives['dD'][0.0428226997],
                            results['dD']['1064nm'],
                            skip_frequency_test=True)

    def test_cook_F_gaussian_b2plyp(self):
        """Check if B2PLYP data are similar to MP2 ones"""

        self.unzip_it(self.zip_F_b2plyp, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_F_b2plyp')
        path = os.path.join(directory, 'nachos_recipe.yml')

        r = files.Recipe(directory=directory)

        with open(path) as f:
            r.read(f)

        fields = preparing.fields_needed_by_recipe(r)

        c = cooking.Cooker(r, directory)
        storage = c.cook([directory])

        # write and read
        self.assertEqual(storage.check(), ([], []))

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
                results = storage.results[t_fields]

                energies = fx.property('computed_energies')

                self.assertNotAlmostEqual(energies['SCF/DFT'], energies['MP2'])
                self.assertAlmostEqual(energies['MP2'], results[''].components[0])

    def test_cook_G_gaussian(self):
        self.unzip_it(self.zip_G, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_G')

        path = os.path.join(directory, 'nachos_recipe.yml')

        r = files.Recipe(directory=directory)

        with open(path) as f:
            r.read(f)

        fields = preparing.fields_needed_by_recipe(r)

        c = cooking.Cooker(r, directory)
        storage = c.cook([directory])

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
                self.assertAlmostEqual(fx.property('computed_energies')['total'], results[''].components[0])

                if level <= 1:
                    electrical_derivatives = fx.property('electrical_derivatives')
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['F']['static'], results['F']['static'])
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['FF']['static'], results['FF']['static'])
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['dD'][0.0428227067], results['dD']['1064nm'], skip_frequency_test=True)

            path = os.path.join(directory, r['name'] + '_{:04d}{}.fchk').format(n, 'b' if level < 2 else '')
            self.assertTrue(os.path.exists(path), msg=path)
            with open(path) as f:
                fx = gaussian.FCHK()
                fx.read(f)

                if level <= 1:
                    geometrical_derivatives = fx.property('geometrical_derivatives')
                    self.assertTensorsAlmostEqual(
                        geometrical_derivatives['G'],
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
        storage = c.cook([directory])
        self.assertEqual(storage.check(), ([], []))

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

                results = storage.results[t_fields]
                self.assertAlmostEqual(fx.property('computed_energies')['total'], results[''].components[0])

                if level <= 1:
                    electrical_derivatives = fx.property('electrical_derivatives')
                    fr = 0.04282271
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['F']['static'], results['F']['static'], skip_frequency_test=True)
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['FF']['static'], results['FF']['static'], skip_frequency_test=True)
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['dD'][fr], results['dD']['1064nm'], skip_frequency_test=True)
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['FFF']['static'], results['FFF']['static'], skip_frequency_test=True)
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['dDF'][fr], results['dDF']['1064nm'], skip_frequency_test=True)
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['XDD'][fr], results['XDD']['1064nm'], skip_frequency_test=True)
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['FFFF']['static'], results['FFFF']['static'], skip_frequency_test=True)
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['dFFD'][fr], results['dFFD']['1064nm'], skip_frequency_test=True)
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['XDDF'][fr], results['XDDF']['1064nm'], skip_frequency_test=True)
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['dDDd'][fr], results['dDDd']['1064nm'], skip_frequency_test=True)
                    self.assertTensorsAlmostEqual(
                        electrical_derivatives['XDDD'][fr], results['XDDD']['1064nm'], skip_frequency_test=True)

    def test_cook_F_qchem(self):
        self.unzip_it(self.zip_F_qchem, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_F_qchem')
        path = os.path.join(directory, 'nachos_recipe.yml')

        r = files.Recipe(directory=directory)

        with open(path) as f:
            r.read(f)

        fields = preparing.fields_needed_by_recipe(r)

        c = cooking.Cooker(r, directory)
        storage = c.cook([directory])
        self.assertEqual(storage.check(), ([], []))

        # check data
        for _ in range(10):
            n = random.randrange(1, len(fields) + 1)
            fields_n, level = fields[n - 1]
            t_fields = tuple(fields_n)
            path = os.path.join(directory, r['name'] + '_{:04d}.log').format(n)
            self.assertTrue(os.path.exists(path), msg=path)

            with open(path) as f:
                fx = nachos.qcip_tools_ext.qchem.QChemLogFile()
                fx.read(f)
                results = storage.results[t_fields]

                self.assertAlmostEqual(fx.property('computed_energies')['total'], results[''].components[0])

    def test_cook_F_scs_mp2(self):
        """Check that using SCS-MP2 is ok"""

        # TODO: relaunch calculation with the inverse E-fields

        self.unzip_it(self.zip_F_scs_mp2, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_F_SCS-MP2')
        path = os.path.join(directory, 'nachos_recipe.yml')

        r = files.Recipe(directory=directory)

        with open(path) as f:
            r.read(f)

        fields = preparing.fields_needed_by_recipe(r)

        c = cooking.Cooker(r, directory)
        storage = c.cook([directory], use_gaussian_logs=True)
        self.assertEqual(storage.check(), ([], []))

        # check data
        for _ in range(10):
            n = random.randrange(1, len(fields) + 1)
            fields_n, level = fields[n - 1]
            t_fields = tuple(-x for x in fields_n)
            path = os.path.join(directory, r['name'] + '_{:04d}.log').format(n)
            self.assertTrue(os.path.exists(path), msg=path)
            with open(path) as f:
                fx = gaussian.Output()
                fx.read(f)

                # check fields
                self.assertEqual(
                    numerical_differentiation.real_fields(fields_n, r['min_field'], r['ratio']),
                    list(fx.property('n:input_electric_field')[1:])
                )

                # check SCS-MP2 energy
                results = storage.results[t_fields]
                energies = fx.property('computed_energies')
                sc_energies = fx.property('n:spin_components_e2')

                self.assertAlmostEqual(
                    results[''].components[0],
                    energies['HF'] + 1 / 3 * (sc_energies['aa'] + sc_energies['bb']) + 6 / 5 * sc_energies['ab'],
                    places=10
                )

    def test_nachos_cook(self):
        """Test the cooker program"""

        self.unzip_it(self.zip_F, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_F')
        recipe_path = os.path.join(directory, 'nachos_recipe.yml')

        h5_path = os.path.join(directory, 'xxx.h5')
        self.assertFalse(os.path.exists(h5_path))

        # process without copy
        process = self.run_python_script(
            'nachos/cook.py', ['-r', recipe_path, '-o', h5_path],
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
