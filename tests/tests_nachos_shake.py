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

    def _test_contributions(self, shaker, test_on, frequencies, are_in, are_not_in=None):
        vibs = shaker.shake(
            frequencies=frequencies,
            only=[(derivatives.Derivative(d[0]), d[1]) for d in test_on])

        self.assertTrue(len(vibs), len(test_on))

        for i, max_level in test_on:
            self.assertIn(i, vibs)

            c = vibs[i]
            self.assertEqual(len(c), len(are_in) + 2)  # the "+2" comes from the two ZPVA contributions

            for j in are_in:
                self.assertIn(j, c)

                if 'D' in i:
                    self.assertEqual(len(c[j]), len(frequencies))
                    for f in frequencies:
                        self.assertIn(f, c[j])
                        self.assertIsInstance(c[j][f], derivatives.Tensor)
                        self.assertEqual(i, c[j][f].representation)
                else:
                    self.assertEqual(len(c[j]), 1)
                    self.assertIn('static', c[j])
                    self.assertIsInstance(c[j]['static'], derivatives.Tensor)
                    self.assertEqual(i, c[j]['static'].representation)

            if are_not_in:
                for j in are_not_in:
                    self.assertNotIn(j, c)

    def test_vibrational_contribution(self):
        """Test the class VibrationalContribution"""

        vcs = [
            # ZPVA
            (('FDF',), 1, 0, ('NNFDF',)),
            (('FDF',), 0, 1, ('NFDF', 'NNN')),
            # pv contrib to beta
            (('F', 'F', 'F'), 1, 0, ('NF', 'NNF')),
            (('F', 'F', 'F'), 0, 1, ('NF', 'NNN')),
            (('F', 'FF'), 0, 0, ('NF', 'NFF')),
            (('F', 'FF'), 1, 1, ('NF', 'NFF', 'NNF', 'NNFF', 'NNN')),
            (('F', 'FF'), 2, 0, ('NNF', 'NNFF')),  # with limitation of the electrical anharmonicity
            (('F', 'FF'), 0, 2, ('NF', 'NFF', 'NNN')),  # with limitation of the electrical anharmonicity
        ]

        for derivs, m, n, needed in vcs:
            vc = shaking.VibrationalContribution(derivs, m, n)
            n = vc.derivatives_needed()

            self.assertEqual(len(needed), len(n), msg=vc)
            for nx in needed:
                self.assertIn(nx, n)

    def test_shaking_alpha_beta(self):
        """Test the shaking class on data from gaussian"""

        df = chemistry_datafile.ChemistryDataFile()

        with open(self.datafile) as f:
            df.read(f)

        shaker = shaking.Shaker(datafile=df)

        # Check everything is there:
        to_check = [
            # ZPVA
            (('F',), 1, 0),
            (('F',), 0, 1),
            (('FD',), 1, 0),
            (('FD',), 0, 1),
            (('FDF',), 1, 0),
            (('FDF',), 0, 1),
            (('FDD',), 1, 0),
            (('FDD',), 0, 1),
            # alpha
            (('F', 'F'), 0, 0),
            (('F', 'F'), 1, 1),
            (('F', 'F'), 2, 0),
            (('F', 'F'), 0, 2),
            # beta:
            (('F', 'F', 'F'), 1, 0),
            (('F', 'F', 'F'), 0, 1),
            (('F', 'FF'), 0, 0),
            (('F', 'FF'), 1, 1),
            (('F', 'FF'), 2, 0),
            (('F', 'FF'), 0, 2),
        ]

        for derivs, m, n in to_check:
            vc = shaking.VibrationalContribution(derivs, m, n)
            self.assertTrue(shaker.check_availability(vc), msg=vc.to_string(fancy=True))

        # test contributions:
        self._test_contributions(
            shaker,
            [('FF', 2), ('FD', 2)],
            [0.02, 0.04],
            ['total', 'total_zpva', 'total_pv', 'F_F__0_0', 'F_F__1_1', 'F_F__2_0', 'F_F__0_2'])

        self._test_contributions(
            shaker,
            [('FF', 1), ('FD', 1)],
            [0.02],
            ['total', 'total_zpva', 'total_pv', 'F_F__0_0'],
            ['F_F__1_1', 'F_F__2_0', 'F_F__0_2'])

        self._test_contributions(
            shaker,
            [('FFF', 2)],
            [],
            ['total', 'total_zpva', 'total_pv', 'F_FF__0_0', 'F_F_F__1_0', 'F_F_F__0_1', 'F_FF__1_1',
             'F_FF__2_0', 'F_FF__0_2'])

        self._test_contributions(
            shaker,
            [('FFF', 1), ('FDF', 1), ('FDD', 1)],
            [0.02],
            ['total', 'total_zpva', 'total_pv', 'F_FF__0_0', 'F_F_F__1_0', 'F_F_F__0_1'],
            ['F_FF__1_1', 'F_FF__2_0', 'F_FF__0_2'])
