from qcip_tools import derivatives, derivatives_e
from qcip_tools.chemistry_files import chemistry_datafile

from nachos.core import analyzing, shaking

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

        # test component access
        component = (1, 1)
        e = analyzing.GetPropertyOfTensor(
            analyzing.component_access,
            explain='whatever',
            component=component)

        self.assertEqual(e.execute(t), t.components[component])

        # dipole does not behave like others!
        t = factories.FakeElectricDipole()
        tx = derivatives.Tensor(t.representation, components=t.components)

        e = analyzing.GetPropertyOfTensor(
            analyzing.property_access,
            analyzing.get_tensor_converter(derivatives_e.ElectricDipole),
            explain='whatever',
            function='norm')

        self.assertEqual(e.execute(tx), t.norm())

    def test_analyzing(self):
        """Test the analysis"""

        df = chemistry_datafile.ChemistryDataFile()

        with open(self.datafile) as f:
            df.read(f)

        vibs = shaking.load_vibrational_contributions(self.datafile, df.spacial_dof)

        analyzer = analyzing.Analyzer(df, vibs)
        analyzer.analyze({
            2: [
                analyzing.GetPropertyOfTensor(
                    analyzing.property_access,
                    analyzing.get_tensor_converter(derivatives_e.PolarisabilityTensor),
                    explain='isotropic_value',
                    function='isotropic_value')
            ],
            3: [
                analyzing.GetPropertyOfTensor(
                    analyzing.component_access,
                    explain='beta_xyz',
                    component=(0, 1, 2))
            ]
        })
