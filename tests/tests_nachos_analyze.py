import subprocess

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

    def test_get_contributions(self):
        self.assertEqual(analyzing.get_vibrational_contributions(2), ['F_F'])
        self.assertEqual(analyzing.get_vibrational_contributions(3), ['F_FF', 'F_F_F'])
        self.assertEqual(analyzing.get_vibrational_contributions(4), ['FF_FF', 'F_FFF', 'F_F_FF', 'F_F_F_F'])

    def test_nachos_analyze(self):
        """Test the analyze command"""

        # base
        process = self.run_python_script(
            'nachos/analyze.py', ['-d', self.datafile, '-O', 'FFF', '-p', 'b::xyz'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        values = [float(a) for a in stdout_t.decode().splitlines()[-2].split()[1:]]

        self.assertAlmostEqual(values[3], sum(values[1:3]), places=4)  # ZPVA
        self.assertAlmostEqual(values[10], sum(values[4:10]), places=4)  # PV
        self.assertAlmostEqual(values[11], sum(values[3:10]), places=4)  # vibs
        self.assertAlmostEqual(values[12], values[11] + values[0], places=4)  # total

        # inverse
        process = self.run_python_script(
            'nachos/analyze.py', ['-d', self.datafile, '-O', 'FFF', '-p', 'b::xyz', '-I'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        values2 = [float(a) for a in stdout_t.decode().splitlines()[-2].split()[1:]]

        self.assertAlmostEqual(values2[0], values[0], places=4)  # elec

        for i in range(1, 10):
            self.assertAlmostEqual(values2[i], values[12] - values[i], places=4, msg=i)

        self.assertAlmostEqual(values2[11], values[11], places=4)  # vib
        self.assertAlmostEqual(values2[12], values[12], places=4)  # tot

        # group
        process = self.run_python_script(
            'nachos/analyze.py', ['-d', self.datafile, '-O', 'FFF', '-p', 'b::xyz', '-g'],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()
        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        values2 = [float(a) for a in stdout_t.decode().splitlines()[-2].split()[1:]]

        self.assertAlmostEqual(values2[0], values[0], places=4)  # elec
        self.assertAlmostEqual(values2[1], values[3], places=4)  # ZPVA
        self.assertAlmostEqual(values2[2], values[4], places=4)  # [µa]⁰
        self.assertAlmostEqual(values2[3], sum(values[5:8]), places=4)  # [µa]ᴵᴵ
        self.assertAlmostEqual(values2[4], sum(values[8:10]), places=4)  # [µ³]ᴵ
        self.assertAlmostEqual(values2[5], values[10], places=4)  # pv
        self.assertAlmostEqual(values2[6], values[11], places=4)  # vib
        self.assertAlmostEqual(values2[7], values[12], places=4)  # tot
