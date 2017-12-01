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
            self.assertTrue(shaker.check_availability(vc), msg=vc)

        # test contributions:
        frequencies = [0.02, 0.04, 0.06]
        p1 = shaker.compute_zpva('10', derivatives.Derivative('FD'), frequencies)
        p2 = shaker.compute_zpva('01', derivatives.Derivative('FD'), frequencies)

        for f in frequencies:
            print('{:<8}'.format(f), '{: .5e} {: .5e}'.format(p1[f].components[0, 0], p2[f].components[0, 0]))

        print('-')
        frequencies = ['static', 0.02, 0.04, 0.06]
        p1 = shaker.compute_pv('mu3_10', derivatives.Derivative('FDF'), frequencies=frequencies)
        p2 = shaker.compute_pv('mu3_10', derivatives.Derivative('FDD'), frequencies=frequencies)

        for f in frequencies:
            print('{:<8}'.format(f), '{: .5e} {: .5e}'.format(p1[f].components[0, 1, 2], p2[f].components[0, 1, 2]))
