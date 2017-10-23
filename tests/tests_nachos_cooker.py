import numpy
import os
import random

from qcip_tools import derivatives, quantities, derivatives_e
from qcip_tools.chemistry_files import gaussian

from tests import NachosTestCase, factories
from nachos.core import files, cooking


class FilesTestCase(NachosTestCase):

    def setUp(self):
        self.geometry = self.copy_to_temporary_directory('water.xyz')
        self.basis_set = self.copy_to_temporary_directory('sto-3g.gbs')
        self.working_directory = self.setup_temporary_directory()

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
        self.assertEqual(len(cooking.fields_needed_by_recipe(r)), 137)
        self.assertEqual(cooking.fields_needed_by_recipe(r)[0], ([0, 0, 0], 1))  # zero field is the first one

        r['k_max'] = 3
        self.assertEqual(len(cooking.fields_needed_by_recipe(r)), 85)

        r['accuracy_level'] = 0
        self.assertEqual(len(cooking.fields_needed_by_recipe(r)), 61)

        r['differentiation'] = {1: ['energy']}
        r._update(r.recipe)
        self.assertEqual(len(cooking.fields_needed_by_recipe(r)), 19)

        r['type'] = 'G'
        r['differentiation'] = {2: ['energy']}
        self.assertEqual(len(cooking.fields_needed_by_recipe(r)), 3 * len(r.geometry) * 2 * r['k_max'] + 1)
        self.assertEqual(cooking.fields_needed_by_recipe(r)[0], ([0] * 3 * len(r.geometry), 1))  # zero field

        r['k_max'] = 5
        r._update(r.recipe)
        self.assertEqual(len(cooking.fields_needed_by_recipe(r)), (3 * len(r.geometry)) ** 2 * 2 * r['k_max'] + 1)

    def test_numerical_differentiation(self):
        energy = 150.
        mu = factories.FakeElectricDipole()
        alpha = factories.FakePolarizabilityTensor(input_fields=(0,))
        beta = factories.FakeFirstHyperpolarizabilityTensor(input_fields=(0, 0))

        def energy_exp(fields, h0, basis, inverse, component, frequency, recipe):
            """Taylor series of the energy"""

            r_field = cooking.Cooker.real_fields(fields, h0, recipe['ratio'])

            x = energy
            x += numpy.tensordot(mu.components, r_field, axes=1)
            x += 1 / 2 * numpy.tensordot(numpy.tensordot(alpha.components, r_field, axes=1), r_field, axes=1)
            x += 1 / 6 * numpy.tensordot(
                numpy.tensordot(numpy.tensordot(beta.components, r_field, axes=1), r_field, axes=1), r_field, axes=1)

            return x

        def dipole_exp(fields, h0, basis, inverse, component, frequency, recipe):
            """Taylor series of the dipole moment"""

            r_field = cooking.Cooker.real_fields(fields, h0, recipe['ratio'])

            x = mu.components.copy()
            x += numpy.tensordot(alpha.components, r_field, axes=1)
            x += 1 / 2 * numpy.tensordot(numpy.tensordot(beta.components, r_field, axes=1), r_field, axes=1)

            return x[component]

        opt_dict = dict(
            flavor='gaussian',
            type='F',
            method='HF',
            basis_set='STO-3G',
            geometry=self.geometry,
            differentiation={3: ['energy', 'F', 'FD']},
            frequencies=['1064nm'],
            accuracy_level=1,
            k_max=5,
            ratio=2,
            min_field=.0004
        )

        r = files.Recipe(**opt_dict)
        r.check_data()

        # compute polarizability
        t, triangles = cooking.compute_numerical_derivative_of_tensor(
            r,
            derivatives.Derivative(from_representation='F', spacial_dof=r.dof),
            derivatives.Derivative('F', spacial_dof=r.dof),
            dipole_exp)

        self.assertArrayAlmostEqual(alpha.components, t.components, places=3)

        t, triangles = cooking.compute_numerical_derivative_of_tensor(
            r,
            derivatives.Derivative(from_representation='', spacial_dof=r.dof),
            derivatives.Derivative('FF', spacial_dof=r.dof),
            energy_exp)

        self.assertArrayAlmostEqual(alpha.components, t.components, places=3)

        # compute first polarizability
        t, triangles = cooking.compute_numerical_derivative_of_tensor(
            r,
            derivatives.Derivative(from_representation='F', spacial_dof=r.dof),
            derivatives.Derivative('FF', spacial_dof=r.dof),
            dipole_exp)

        self.assertArrayAlmostEqual(beta.components, t.components, delta=.001)

        t, triangles = cooking.compute_numerical_derivative_of_tensor(
            r,
            derivatives.Derivative(from_representation='', spacial_dof=r.dof),
            derivatives.Derivative('FFF', spacial_dof=r.dof),
            energy_exp)

        self.assertArrayAlmostEqual(beta.components, t.components, delta=.01)

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
        deformed = cooking.Cooker.deform_geometry(r.geometry, zero_fields)
        for i, a in enumerate(deformed):
            self.assertArrayAlmostEqual(a.position, r.geometry[i].position)

        # simple deformation:
        fields = zero_fields.copy()
        fields[0] = 1  # x of first atom
        deformed = cooking.Cooker.deform_geometry(
            r.geometry, cooking.Cooker.real_fields(fields, r['min_field'], r['ratio']))

        self.assertNotEqual(r.geometry[0].position[0], deformed[0].position[0])
        self.assertEqual(r.geometry[0].position[1], deformed[0].position[1])  # only the first coordinate was modified
        self.assertEqual(r.geometry[0].position[2], deformed[0].position[2])

        for i, a in enumerate(deformed[1:]):
            self.assertArrayAlmostEqual(a.position, r.geometry[i + 1].position)  # nothing else was deformed !

        self.assertEqual(r.geometry[0].position[0] + r['min_field'] * quantities.AuToAngstrom, deformed[0].position[0])

        # double deformation:
        fields = zero_fields.copy()
        fields[0] = 3
        fields[4] = -2  # y of second atom
        deformed = cooking.Cooker.deform_geometry(
            r.geometry, cooking.Cooker.real_fields(fields, r['min_field'], r['ratio']))

        self.assertNotEqual(r.geometry[0].position[0], deformed[0].position[0])
        self.assertNotEqual(r.geometry[1].position[1], deformed[1].position[1])

        self.assertEqual(
            r.geometry[0].position[0] + r['min_field'] * r['ratio'] ** 2 * quantities.AuToAngstrom,
            deformed[0].position[0])

        self.assertEqual(
            r.geometry[1].position[1] - r['min_field'] * r['ratio'] * quantities.AuToAngstrom,
            deformed[1].position[1])

    def test_cooker_for_F(self):
        """Test the cooker class"""

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

        fields = cooking.fields_needed_by_recipe(r)

        cook = cooking.Cooker(recipe=r, directory=self.working_directory)
        cook.cook()

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
                    cooking.Cooker.real_fields(fields_n, min_field, 2.))

    def test_cooker_for_G(self):
        """Test the cooker class"""

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

        fields = cooking.fields_needed_by_recipe(r)

        recipe_path = os.path.join(self.working_directory, 'nachos_recipe.yaml')
        with open(recipe_path, 'w') as f:
            r.write(f)

        cook = cooking.Cooker(recipe=r, directory=self.working_directory)
        cook.cook()

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
                self.assertArrayAlmostEqual(r.geometry[i].position, a.position)  # geometry not modified, it is base

        path = os.path.join(self.working_directory, name + '_0001b.com')
        self.assertTrue(os.path.exists(path))
        with open(path) as f:
            fi = gaussian.Input()
            fi.read(f)

            self.assertFalse(any(['polar' in a for a in fi.input_card]))
            self.assertEqual(len(fi.other_blocks), 1)
            self.assertEqual(fi.other_blocks[0][0], 'O     0')

            for i, a in enumerate(fi.molecule):
                self.assertArrayAlmostEqual(r.geometry[i].position, a.position)  # geometry not modified, it is base

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

                deformed = cooking.Cooker.deform_geometry(
                    r.geometry, cooking.Cooker.real_fields(fields_n, min_field, 2))

                for i, a in enumerate(fi.molecule):
                    self.assertArrayAlmostEqual(deformed[i].position, a.position)

    def test_nachos_cooker(self):
        """Test the cooker program"""
        pass
