from qcip_tools import derivatives, derivatives_e

from nachos.core import analyzing

from tests import NachosTestCase, factories


class AnalyzeTestCase(NachosTestCase):

    def setUp(self):
        self.datafile = self.copy_to_temporary_directory('molecule_with_derivs_and_vibs.h5')

    def test_property_get(self):
        """Test if one can use property accessor"""

        t = factories.FakePolarizabilityTensor()
        tx = derivatives.Tensor(
            t.representation, frequency=t.frequency, spacial_dof=t.spacial_dof, components=t.components)

        # test property access
        e = analyzing.GetPropertyOfTensor(
            analyzing.property_access,
            analyzing.get_tensor_converter(derivatives_e.PolarisabilityTensor),
            explain='whatever',
            function='isotropic_value')

        self.assertEqual(e.execute(tx), t.isotropic_value())

        # test inverse
        self.assertEqual(e.execute(tx, total=tx), 0.0)

        # test component access
        component = (1, 1)
        e = analyzing.GetPropertyOfTensor(
            analyzing.component_access,
            explain='whatever',
            component=component)

        self.assertEqual(e.execute(t), t.components[component])

        # test inverse
        self.assertEqual(e.execute(tx, total=tx), 0.0)

        # dipole does not behave like others!
        t = factories.FakeElectricDipole()
        tx = derivatives.Tensor(t.representation, components=t.components)

        e = analyzing.GetPropertyOfTensor(
            analyzing.property_access,
            analyzing.get_tensor_converter(derivatives_e.ElectricDipole),
            explain='whatever',
            function='norm')

        self.assertEqual(e.execute(tx), t.norm())
