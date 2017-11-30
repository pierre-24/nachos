from qcip_tools import derivatives
from qcip_tools.chemistry_files import chemistry_datafile

from nachos.core import shaking

from tests import NachosTestCase


class ShakeTestCase(NachosTestCase):

    def setUp(self):
        self.datafile = self.copy_to_temporary_directory('molecule_with_derivs.h5')

    def tearDown(self):
        super().tearDown()
        pass

    def test_shaking_alpha_beta(self):
        """Test the shaking class on data from gaussian"""

        df = chemistry_datafile.ChemistryDataFile()

        with open(self.datafile) as f:
            df.read(f)

        shaker = shaking.Shaker(datafile=df)

        # Check everything is there:
        to_ZPVA = ['F', 'FF', 'FD', 'FFF', 'FDF', 'FDD']
        for t in to_ZPVA:
            self.assertTrue(shaker.check_availability((derivatives.Derivative(t),), True, 1, 0), msg=t)  # []¹⁰
            self.assertTrue(shaker.check_availability((derivatives.Derivative(t),), True, 0, 1), msg=t)  # []⁰¹

        self.assertTrue(shaker.check_availability((derivatives.Derivative('F'),), False, 0, 0))  # [µ²]⁰⁰
        self.assertTrue(shaker.check_availability((derivatives.Derivative('F'),), False, 1, 1))  # [µ²]¹¹
        self.assertTrue(shaker.check_availability((derivatives.Derivative('F'),), False, 2, 0))  # [µ²]²⁰
        self.assertTrue(shaker.check_availability((derivatives.Derivative('F'),), False, 0, 2))  # [µ²]⁰²

        self.assertTrue(shaker.check_availability(
            (derivatives.Derivative('F'), derivatives.Derivative('FF')), False, 0, 0))  # [µa]⁰⁰
        self.assertTrue(shaker.check_availability(
            (derivatives.Derivative('F'), derivatives.Derivative('FF')), False, 1, 1))  # [µa]¹¹
        self.assertTrue(shaker.check_availability(
            (derivatives.Derivative('F'), derivatives.Derivative('FF')), False, 2, 0))  # [µa]²⁰
        self.assertTrue(shaker.check_availability(
            (derivatives.Derivative('F'), derivatives.Derivative('FF')), False, 0, 2))  # [µa]⁰²

        self.assertTrue(shaker.check_availability((derivatives.Derivative('F'),), False, 1, 0))  # [µ³]¹⁰
        self.assertTrue(shaker.check_availability((derivatives.Derivative('F'),), False, 0, 1))  # [µ³]⁰¹

        # test!
        frequencies = ['static', 0.02, 0.04, 0.06]
        # polarizabilities = shaker.compute_zpva('01', derivatives.Derivative('FD'), frequencies)

        polarizabilities = shaker.compute_pv('mu2_02', derivatives.Derivative('FD'), frequencies=frequencies)

        for f in frequencies:
            print('{:<8}'.format(f), '{: .5e}'.format(polarizabilities[f].isotropic_value()))
