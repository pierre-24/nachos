import os
import math
import subprocess

import numpy

from tests import NachosTestCase

from qcip_tools import derivatives, derivatives_g
from qcip_tools.chemistry_files import gaussian, chemistry_datafile

from nachos.core import files, baking


class BakeTestCase(NachosTestCase):

    def setUp(self):
        self.zip_F = 'numdiff_F.zip'
        self.zip_G = 'numdiff_G.zip'
        self.zip_G_dalton = 'numdiff_G_dalton.zip'
        self.working_directory = self.setup_temporary_directory()

    def tearDown(self):
        super().tearDown()
        pass

    def test_bake_gaussian_F(self):
        self.unzip_it(self.zip_F, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_F')
        recipe_path = os.path.join(directory, 'nachos_recipe.yml')
        storage_path = os.path.join(directory, 'verification', 'nachos_data.h5')
        solution_path = os.path.join(directory, 'verification', 'water_test_0001.fchk')

        r = files.Recipe(directory=directory)

        with open(recipe_path) as f:
            r.read(f)

        fchk = gaussian.FCHK()
        with open(solution_path) as f:
            fchk.read(f)

        electrical_derivatives = fchk.property('electrical_derivatives')

        storage = files.ComputationalResults(r, directory=directory)
        storage.read(storage_path)
        baker = baking.Baker(r, storage, directory=directory)

        # with energy:
        cf_with_energy = baker.bake(only=[(derivatives.Derivative(), 3)])

        self.assertIn('F', cf_with_energy.derivatives)
        self.assertIn('FF', cf_with_energy.derivatives)
        self.assertIn('FFF', cf_with_energy.derivatives)
        self.assertEqual(len(cf_with_energy.derivatives), 3)

        self.assertTensorsAlmostEqual(
            electrical_derivatives['F']['static'], cf_with_energy.derivatives['F']['static'])

        self.assertTensorsAlmostEqual(
            electrical_derivatives['FF']['static'], cf_with_energy.derivatives['FF']['static'])

        self.assertTensorsAlmostEqual(
            electrical_derivatives['FFF']['static'], cf_with_energy.derivatives['FFF']['static'])

        # with mu:
        cf_with_mu = baker.bake(only=[(derivatives.Derivative('F'), 2)])

        self.assertIn('FF', cf_with_mu.derivatives)
        self.assertIn('FFF', cf_with_mu.derivatives)
        self.assertEqual(len(cf_with_mu.derivatives), 2)

        self.assertTensorsAlmostEqual(
            electrical_derivatives['FF']['static'], cf_with_mu.derivatives['FF']['static'])

        self.assertTensorsAlmostEqual(
            electrical_derivatives['FFF']['static'],
            cf_with_mu.derivatives['FFF']['static'],
            delta=.01)

        # with alpha:
        cf_with_alpha = baker.bake(only=[(derivatives.Derivative('FF'), 1)])

        self.assertIn('FFF', cf_with_alpha.derivatives)
        self.assertEqual(len(cf_with_alpha.derivatives), 1)

        self.assertTensorsAlmostEqual(
            electrical_derivatives['FFF']['static'],
            cf_with_alpha.derivatives['FFF']['static'])

        # dynamic
        cf_with_alpha = baker.bake(only=[(derivatives.Derivative('dD'), 1)])

        self.assertIn('dDF', cf_with_alpha.derivatives)
        self.assertEqual(len(cf_with_alpha.derivatives), 1)

        self.assertTensorsAlmostEqual(
            electrical_derivatives['dDF'][0.0428226997],
            cf_with_alpha.derivatives['dDF']['1064nm'],
            skip_frequency_test=True)

        # fire some errors:
        with self.assertRaises(baking.BadBaking):
            baker.bake(only=[(derivatives.Derivative(), 4)])

        with self.assertRaises(baking.BadBaking):
            baker.bake(only=[(derivatives.Derivative('dDF'), 0)])

    def test_bake_gaussian_F_with_force(self):
        self.unzip_it(self.zip_F, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_F')
        recipe_path = os.path.join(directory, 'nachos_recipe.yml')
        storage_path = os.path.join(directory, 'verification', 'nachos_data.h5')
        solution_path = os.path.join(directory, 'verification', 'water_test_0001.fchk')

        r = files.Recipe(directory=directory)

        with open(recipe_path) as f:
            r.read(f)

        fchk = gaussian.FCHK()
        with open(solution_path) as f:
            fchk.read(f)

        electrical_derivatives = fchk.property('electrical_derivatives')

        storage = files.ComputationalResults(r, directory=directory)
        storage.read(storage_path)
        baker = baking.Baker(r, storage, directory=directory)

        cf_free = baker.bake(only=[(derivatives.Derivative(), 2)])
        cf_force = baker.bake(only=[(derivatives.Derivative(), 2)], force_choice=(2, 2))
        self.assertIn('FF', cf_force.derivatives)

        self.assertTensorsAlmostEqual(
            electrical_derivatives['FF']['static'], cf_force.derivatives['FF']['static'])

        diffs = cf_free.derivatives['FF']['static'].components - cf_force.derivatives['FF']['static'].components
        self.assertTrue(numpy.all(diffs < 1e-3))

    def test_bake_gaussian_G(self):
        self.unzip_it(self.zip_G, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_G')
        recipe_path = os.path.join(directory, 'nachos_recipe.yml')
        storage_path = os.path.join(directory, 'verification', 'nachos_data.h5')
        solution_path = os.path.join(directory, 'verification', 'water_test_0001b.fchk')

        r = files.Recipe(directory=directory)

        with open(recipe_path) as f:
            r.read(f)

        storage = files.ComputationalResults(r, directory=directory)
        storage.read(storage_path)

        fchk = gaussian.FCHK()
        with open(solution_path) as f:
            fchk.read(f)

        geometrical_derivatives = fchk.property('geometrical_derivatives')

        baker = baking.Baker(r, storage, directory=directory)

        # with energy:
        cf_with_energy = baker.bake(only=[(derivatives.Derivative(), 2)])

        self.assertIn('G', cf_with_energy.derivatives)
        self.assertIn('GG', cf_with_energy.derivatives)
        self.assertEqual(len(cf_with_energy.derivatives), 2)

        self.assertTensorsAlmostEqual(
            geometrical_derivatives['G'], cf_with_energy.derivatives['G'])

        self.assertTensorsAlmostEqual(
            geometrical_derivatives['GG'], cf_with_energy.derivatives['GG'])

        # try with all
        cf = baker.bake()

        self.assertIn('G', cf.derivatives)
        self.assertIn('GG', cf.derivatives)
        self.assertIn('GF', cf.derivatives)
        self.assertIn('GFF', cf.derivatives)
        self.assertIn('GdD', cf.derivatives)

        self.assertEqual(len(cf.derivatives), 5)

        self.assertTensorsAlmostEqual(
            geometrical_derivatives['GG'], cf.derivatives['GG'])

        self.assertArraysAlmostEqual(
            cf.derivatives['GF']['static'].components.flatten(), fchk.get('Dipole Derivatives'))

        dalphadq = fchk.get('Derivative Alpha(-w,w)')
        self.assertArraysAlmostEqual(cf.derivatives['GFF']['static'].components.flatten(), dalphadq[:81])
        self.assertArraysAlmostEqual(cf.derivatives['GdD']['1064nm'].components.flatten(), dalphadq[81:])

        # try projection
        mwh = derivatives_g.MassWeightedHessian(fchk.molecule, geometrical_derivatives['GG'].components)
        baking.project_geometrical_derivatives(r, cf, mwh)

        self.assertIn('N', cf.derivatives)
        self.assertIn('NN', cf.derivatives)
        self.assertIn('NF', cf.derivatives)
        self.assertIn('NFF', cf.derivatives)
        self.assertIn('NdD', cf.derivatives)

        ph = cf.derivatives['NN'].components
        for i in range(r.dof):
            # so the square of the frequencies are in the diagonal:
            self.assertAlmostEqual(
                math.fabs(ph[i, i]), mwh.frequencies[i] ** 2, places=5)

        # now, bake and steal results from zero field
        cf_with_copy = baker.bake(only=[(derivatives.Derivative(), 0)], copy_zero_field_basis=True)

        # gradient and hessian ...
        self.assertIn('G', cf_with_copy.derivatives)
        self.assertIn('GG', cf_with_copy.derivatives)

        # ... but also:
        self.assertIn('', cf_with_copy.derivatives)
        self.assertIn('F', cf_with_copy.derivatives)
        self.assertIn('FF', cf_with_copy.derivatives)
        self.assertIn('dD', cf_with_copy.derivatives)

        self.assertEqual(len(cf_with_copy.derivatives), 6)

    def test_bake_dalton_G(self):
        self.unzip_it(self.zip_G_dalton, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_G_dalton')
        recipe_path = os.path.join(directory, 'nachos_recipe.yml')
        storage_path = os.path.join(directory, 'verification', 'nachos_data.h5')

        r = files.Recipe(directory=directory)

        with open(recipe_path) as f:
            r.read(f)

        storage = files.ComputationalResults(r, directory=directory)
        storage.read(storage_path)

        baker = baking.Baker(r, storage, directory=directory)
        cf = baker.bake()

        # only test what was computed, since there is no way to test
        # (G, GG and GF are the same as with HF, but anyway)
        should_be_in = [
            'G', 'GG', 'GF', 'GFF', 'GdD', 'GFFF', 'GdDF', 'GXDD', 'GFFFF', 'GdFFD', 'GXDDF', 'GdDDd', 'GXDDD']

        for i in should_be_in:
            self.assertIn(i, cf.derivatives, msg=i)

        self.assertEqual(len(should_be_in), len(cf.derivatives))

    def test_nachos_bake(self):
        """Test the baking program"""

        self.unzip_it(self.zip_G, self.working_directory)
        directory = os.path.join(self.working_directory, 'numdiff_G')
        recipe_path = os.path.join(directory, 'nachos_recipe.yml')
        storage_path = os.path.join(directory, 'verification', 'nachos_data.h5')

        h5_path = os.path.join(directory, 'out.h5')
        self.assertFalse(os.path.exists(h5_path))

        # process without projection
        process = self.run_python_script(
            'nachos/bake.py', ['-r', recipe_path, '-d', storage_path, '-o', h5_path],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertEqual(len(stdout_t), 0, msg=stdout_t.decode())

        self.assertTrue(os.path.exists(h5_path))

        # check h5 file
        cf = chemistry_datafile.ChemistryDataFile()
        with open(h5_path) as f:
            cf.read(f)

        must_be_in = ['G', 'GG', '', 'F', 'FF', 'dD', 'G', 'GG', 'GF', 'GFF', 'GdD']
        for m in must_be_in:
            self.assertIn(m, cf.derivatives)

        # process with projection
        process = self.run_python_script(
            'nachos/bake.py', ['-r', recipe_path, '-d', storage_path, '-o', h5_path, '-p'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertEqual(len(stdout_t), 0, msg=stdout_t.decode())

        self.assertTrue(os.path.exists(h5_path))

        # check h5 file
        cf = chemistry_datafile.ChemistryDataFile()
        with open(h5_path) as f:
            cf.read(f)

        must_also_be_in = ['N', 'NN', 'NF', 'NFF', 'NdD']
        for m in must_be_in:
            self.assertIn(m, cf.derivatives)
        for m in must_also_be_in:
            self.assertIn(m, cf.derivatives)

        # process with only (and "do-not-steal")
        process = self.run_python_script(
            'nachos/bake.py', ['-r', recipe_path, '-d', storage_path, '-o', h5_path, '-O', 'F;energy:1', '-S'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertEqual(len(stdout_t), 0, msg=stdout_t.decode())

        self.assertTrue(os.path.exists(h5_path))

        # check h5 file
        cf = chemistry_datafile.ChemistryDataFile()
        with open(h5_path) as f:
            cf.read(f)

        must_be_in_x = ['GF', 'G']
        for m in must_be_in_x:
            self.assertIn(m, cf.derivatives)
        self.assertEqual(len(cf.derivatives), 2)
