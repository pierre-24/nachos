import numpy

from qcip_tools import derivatives, quantities

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
            type='electric',
            method='HF',
            basis_set='STO-3G',
            geometry=self.geometry,
            bases=['energy', 'F', 'FD'],
            frequencies=['1064nm'],
            accuracy_level=1,
            k_max=5,
            max_differentiation=3,
            ratio=2,
            min_field=.001
        )

        r = files.Recipe(**opt_dict)
        self.assertEqual(len(cooking.fields_needed_by_recipe(r)), 137)
        self.assertEqual(cooking.fields_needed_by_recipe(r)[0], [0, 0, 0])

        r['k_max'] = 3
        self.assertEqual(len(cooking.fields_needed_by_recipe(r)), 85)

        r['accuracy_level'] = 0
        self.assertEqual(len(cooking.fields_needed_by_recipe(r)), 61)

        r['max_differentiation'] = 1
        self.assertEqual(len(cooking.fields_needed_by_recipe(r)), 19)

        r['type'] = 'geometric'
        self.assertEqual(len(cooking.fields_needed_by_recipe(r)), 3 * len(r.geometry) * 2 * r['k_max'] + 1)
        self.assertEqual(cooking.fields_needed_by_recipe(r)[0], [0] * 3 * len(r.geometry))

        r['max_differentiation'] = 2
        r['k_max'] = 5
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
            type='electric',
            method='HF',
            basis_set='STO-3G',
            geometry=self.geometry,
            bases=['energy', 'F', 'FD'],
            frequencies=['1064nm'],
            accuracy_level=1,
            k_max=5,
            max_differentiation=3,
            ratio=2,
            min_field=.0004
        )

        r = files.Recipe(**opt_dict)
        dof = 3 * len(r.geometry)

        # compute polarizability
        t, triangles = cooking.compute_numerical_derivative_of_tensor(
            r,
            derivatives.Derivative(from_representation='F', spacial_dof=dof),
            derivatives.Derivative('F', spacial_dof=dof),
            dipole_exp)

        self.assertArrayAlmostEqual(alpha.components, t.components, places=3)

        t, triangles = cooking.compute_numerical_derivative_of_tensor(
            r,
            derivatives.Derivative(from_representation='', spacial_dof=dof),
            derivatives.Derivative('FF', spacial_dof=dof),
            energy_exp)

        self.assertArrayAlmostEqual(alpha.components, t.components, places=3)

        # compute first polarizability
        t, triangles = cooking.compute_numerical_derivative_of_tensor(
            r,
            derivatives.Derivative(from_representation='F', spacial_dof=dof),
            derivatives.Derivative('FF', spacial_dof=dof),
            dipole_exp)

        self.assertArrayAlmostEqual(beta.components, t.components, delta=.001)

        t, triangles = cooking.compute_numerical_derivative_of_tensor(
            r,
            derivatives.Derivative(from_representation='', spacial_dof=dof),
            derivatives.Derivative('FFF', spacial_dof=dof),
            energy_exp)

        self.assertArrayAlmostEqual(beta.components, t.components, delta=.01)

    def test_deform_geometry(self):
        """Test geometry deformation"""

        opt_dict = dict(
            flavor='gaussian',
            type='geometric',
            method='HF',
            basis_set='gen',
            geometry=self.geometry,
            bases=['energy', 'F', 'FD'],
            frequencies=['1064nm'],
            min_field=0.01  # atomic units
        )

        r = files.Recipe(**opt_dict)
        r.check_data()
        dof = 3 * len(r.geometry)

        zero_fields = [0] * dof

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
            self.assertArrayAlmostEqual(a.position, r.geometry[i+1].position)  # nothing else was deformed !

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

    def test_cooker(self):
        """Test the cooker class"""

        opt_dict = dict(
            flavor='gaussian',
            type='electric',
            method='HF',
            basis_set='gen',
            geometry=self.geometry,
            bases=['energy', 'F', 'FD'],
            frequencies=['1064nm']
        )

        flavor_extra = {
            'memory': '3Gb',
            'procs': 4,
            'gen_basis': self.basis_set
        }

        r = files.Recipe(flavor_extra=flavor_extra, **opt_dict)
        r.check_data()

        cook = cooking.Cooker(recipe=r, directory=self.working_directory)
        cook.cook()

        print(self.working_directory)

    def test_nachos_cooker(self):
        """Test the cooker program"""
        pass
