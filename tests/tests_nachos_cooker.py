import numpy

from qcip_tools import derivatives, numerical_differentiation

from tests import NachosTestCase, factories
from nachos.core import files, cooking


class FilesTestCase(NachosTestCase):

    def setUp(self):
        self.geometry = self.copy_to_temporary_directory('water.xyz')
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

            r_field = [h0 * numerical_differentiation.ak_shifted(recipe['ratio'], f) for f in fields]

            x = energy
            x += numpy.tensordot(mu.components, r_field, axes=1)
            x += 1 / 2 * numpy.tensordot(numpy.tensordot(alpha.components, r_field, axes=1), r_field, axes=1)
            x += 1 / 6 * numpy.tensordot(
                numpy.tensordot(numpy.tensordot(beta.components, r_field, axes=1), r_field, axes=1), r_field, axes=1)

            return x

        def dipole_exp(fields, h0, basis, inverse, component, frequency, recipe):
            """Taylor series of the dipole moment"""

            r_field = [h0 * numerical_differentiation.ak_shifted(recipe['ratio'], f) for f in fields]

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
            min_field=.0001
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

    def test_cooker(self):
        """Test the cooker class"""

        opt_dict = dict(
            flavor='gaussian',
            type='electric',
            method='HF',
            basis_set='STO-3G',
            geometry=self.geometry,
            bases=['energy', 'F', 'FD'],
            frequencies=['1064nm']
        )

        flavor_extra = {
            'memory': '3Gb',
            'procs': 4
        }

        r = files.Recipe(flavor_extra=flavor_extra, **opt_dict)
        r.check_data()

        cook = cooking.Cooker(recipe=r, directory=self.working_directory)
        cook.cook()

    def test_nachos_cooker(self):
        """Test the cooker program"""
        pass
