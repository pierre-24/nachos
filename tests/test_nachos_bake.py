import os

from tests import NachosTestCase

from qcip_tools import derivatives
from qcip_tools.chemistry_files import gaussian

from nachos.core import files, baking


class BakeTestCase(NachosTestCase):

    def setUp(self):
        self.zip_F = self.copy_to_temporary_directory('numdiff_F.zip')
        self.zip_G = self.copy_to_temporary_directory('numdiff_G.zip')
        self.zip_G_dalton = self.copy_to_temporary_directory('numdiff_G_dalton.zip')
        self.working_directory = self.setup_temporary_directory()

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

        self.assertArrayAlmostEqual(
            electrical_derivatives['F']['static'].components, cf_with_energy.derivatives['F']['static'])

        self.assertArrayAlmostEqual(
            electrical_derivatives['FF']['static'].components, cf_with_energy.derivatives['FF']['static'])

        self.assertArrayAlmostEqual(
            electrical_derivatives['FFF']['static'].components, cf_with_energy.derivatives['FFF']['static'])

        # with mu:
        cf_with_mu = baker.bake(only=[(derivatives.Derivative('F'), 2)])

        self.assertIn('FF', cf_with_mu.derivatives)
        self.assertIn('FFF', cf_with_mu.derivatives)
        self.assertEqual(len(cf_with_mu.derivatives), 2)

        self.assertArrayAlmostEqual(
            electrical_derivatives['FF']['static'].components, cf_with_mu.derivatives['FF']['static'])

        self.assertArrayAlmostEqual(
            electrical_derivatives['FFF']['static'].components,
            cf_with_mu.derivatives['FFF']['static'],
            delta=.01)

        # with alpha:
        cf_with_alpha = baker.bake(only=[(derivatives.Derivative('FF'), 1)])

        self.assertIn('FFF', cf_with_alpha.derivatives)
        self.assertEqual(len(cf_with_alpha.derivatives), 1)

        self.assertArrayAlmostEqual(
            electrical_derivatives['FFF']['static'].components,
            cf_with_alpha.derivatives['FFF']['static'])

        # fire some errors:
        with self.assertRaises(baking.BadBaking):
            baker.bake(only=[(derivatives.Derivative(), 4)])

        with self.assertRaises(baking.BadBaking):
            baker.bake(only=[(derivatives.Derivative('FDF'), 0)])

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

        self.assertArrayAlmostEqual(
            geometrical_derivatives['G'].components, cf_with_energy.derivatives['G'])

        self.assertArrayAlmostEqual(
            geometrical_derivatives['GG'].components, cf_with_energy.derivatives['GG'])

        # try with all
        cf = baker.bake()

        self.assertIn('G', cf.derivatives)
        self.assertIn('GG', cf.derivatives)
        self.assertIn('GF', cf.derivatives)
        self.assertIn('GFF', cf.derivatives)
        self.assertIn('GFD', cf.derivatives)

        self.assertEqual(len(cf.derivatives), 5)

        self.assertArrayAlmostEqual(
            geometrical_derivatives['GG'].components, cf.derivatives['GG'])

        self.assertArrayAlmostEqual(
            cf.derivatives['GF']['static'].flatten(), fchk.get('Dipole Derivatives'))

        dalphadq = fchk.get('Derivative Alpha(-w,w)')
        self.assertArrayAlmostEqual(cf.derivatives['GFF']['static'].flatten(), dalphadq[:81])
        self.assertArrayAlmostEqual(cf.derivatives['GFD']['1064nm'].flatten(), dalphadq[81:])

        # now, bake and steal results from zero field
        cf_with_copy = baker.bake(only=[(derivatives.Derivative(), 0)], copy_zero_field_basis=True)

        # gradient and hessian ...
        self.assertIn('G', cf_with_copy.derivatives)
        self.assertIn('GG', cf_with_copy.derivatives)

        # ... but also:
        self.assertIn('', cf_with_copy.derivatives)
        self.assertIn('F', cf_with_copy.derivatives)
        self.assertIn('FF', cf_with_copy.derivatives)
        self.assertIn('FD', cf_with_copy.derivatives)

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
            'G', 'GG', 'GF', 'GFF', 'GFD', 'GFFF', 'GFDF', 'GFDD', 'GFFFF', 'GFDFF', 'GFDDF', 'GFDDd', 'GFFF']

        for i in should_be_in:
            self.assertIn(i, cf.derivatives, msg=i)

        self.assertEqual(len(should_be_in), len(cf.derivatives))
