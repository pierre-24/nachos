import subprocess

from qcip_tools import derivatives, derivatives_e
from qcip_tools.chemistry_files import chemistry_datafile

from nachos.core import shaking

from tests import NachosTestCase


class ShakeTestCase(NachosTestCase):

    def setUp(self):
        self.datafile = self.copy_to_temporary_directory('molecule_with_derivs.h5')
        self.datafile_g = self.copy_to_temporary_directory('molecule_with_derivs_gamma.h5')

    def tearDown(self):
        super().tearDown()
        pass

    def _test_contributions(self, shaker, test_on, frequencies, are_in, are_not_in=None, is_zpva=True):
        vibs = shaker.shake(
            frequencies=frequencies,
            only=[(derivatives.Derivative(d[0]), d[1]) for d in test_on])

        self.assertTrue(len(vibs), len(test_on))

        for i, max_level in test_on:
            self.assertIn(i, vibs)

            c = vibs[i]
            self.assertEqual(len(c.vibrational_contributions), len(are_in) + (2 if is_zpva else 0))
            # the "+2" comes from the two ZPVA contributions

            for j in are_in:
                self.assertIn(j, c.vibrational_contributions)
                vc = shaking.VibrationalContribution.from_representation(j)
                if vc.zpva:
                    self.assertIn(vc, c.per_type['zpva'])
                else:
                    self.assertIn(vc, c.per_type['pv'])

                if 'D' in i:
                    self.assertEqual(len(c.vibrational_contributions[j]), len(frequencies))
                    for f in frequencies:
                        if vc.zpva:
                            self.assertIn(f, c.frequencies_zpva)
                        else:
                            self.assertIn(f, c.frequencies_pv)

                        self.assertIn(f, c.vibrational_contributions[j])
                        self.assertIsInstance(c.vibrational_contributions[j][f], derivatives.Tensor)
                        self.assertEqual(i, c.vibrational_contributions[j][f].representation)
                else:
                    self.assertEqual(len(c.vibrational_contributions[j]), 1)
                    self.assertIn('static', c.vibrational_contributions[j])
                    self.assertIsInstance(c.vibrational_contributions[j]['static'], derivatives.Tensor)
                    self.assertEqual(i, c.vibrational_contributions[j]['static'].representation)

            if are_not_in:
                for j in are_not_in:
                    self.assertNotIn(j, c.vibrational_contributions)

        return vibs

    def test_vibrational_contribution(self):
        """Test the class VibrationalContribution"""

        vcs = [
            # ZPVA
            (('dDF',), 1, 0, ('NNdDF',)),
            (('dDF',), 0, 1, ('NdDF', 'NNN')),
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

            # test from representation
            vcx = shaking.VibrationalContribution.from_representation(vc.to_string())
            self.assertEqual(vcx.derivatives, vc.derivatives)
            self.assertEqual((vcx.m, vcx.n), (vc.m, vc.n))

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
            (('dD',), 1, 0),
            (('dD',), 0, 1),
            (('dDF',), 1, 0),
            (('dDF',), 0, 1),
            (('XDD',), 1, 0),
            (('XDD',), 0, 1),
            # alpha pv
            (('F', 'F'), 0, 0),
            (('F', 'F'), 1, 1),
            (('F', 'F'), 2, 0),
            (('F', 'F'), 0, 2),
            # beta pv
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
            [('FF', 2), ('dD', 2)],
            [0.02, 0.04],
            ['F_F__0_0', 'F_F__1_1', 'F_F__2_0', 'F_F__0_2'])

        self._test_contributions(
            shaker,
            [('FF', 1), ('dD', 1)],
            [0.02],
            ['F_F__0_0'],
            ['F_F__1_1', 'F_F__2_0', 'F_F__0_2'])

        self._test_contributions(
            shaker,
            [('FFF', 2)],
            [],
            ['F_FF__0_0', 'F_F_F__1_0', 'F_F_F__0_1', 'F_FF__1_1',
             'F_FF__2_0', 'F_FF__0_2'])

        self._test_contributions(
            shaker,
            [('FFF', 1), ('dDF', 1), ('XDD', 1)],
            [0.02],
            ['F_FF__0_0', 'F_F_F__1_0', 'F_F_F__0_1'],
            ['F_FF__1_1', 'F_FF__2_0', 'F_FF__0_2'])

        self._test_contributions(
            shaker,
            [('FFFF', 1), ('dDFF', 1), ('XDDF', 1)],
            [0.02],
            ['FF_FF__0_0', 'F_FFF__0_0', 'F_F_FF__1_0', 'F_F_FF__0_1'], is_zpva=False)

    def test_shaking_gamma(self):
        """Test the shaking class on data from a dalton that contains gamma"""

        df = chemistry_datafile.ChemistryDataFile()

        with open(self.datafile_g) as f:
            df.read(f)

        shaker = shaking.Shaker(datafile=df)

        # Check everything is there:
        to_check = [
            # ZPVA
            (('FFFF',), 1, 0),
            (('FFFF',), 0, 1),
            (('XDDD',), 1, 0),
            (('XDDD',), 0, 1),
            # gamma pv
            (('F', 'F', 'FF'), 1, 0),
            (('F', 'F', 'FF'), 0, 1),
            (('FF', 'FF'), 0, 0),
            (('FF', 'FF'), 1, 1),
            (('FF', 'FF'), 2, 0),
            (('FF', 'FF'), 0, 2),
            (('F', 'FFF'), 0, 0),
            (('F', 'FFF'), 1, 1),
            (('F', 'FFF'), 2, 0),
            (('F', 'FFF'), 0, 2),
            (('F', 'F', 'F', 'F'), 1, 1),
            (('F', 'F', 'F', 'F'), 2, 0),
            (('F', 'F', 'F', 'F'), 0, 2),
        ]

        for derivs, m, n in to_check:
            vc = shaking.VibrationalContribution(derivs, m, n)
            self.assertTrue(shaker.check_availability(vc), msg=vc.to_string(fancy=True))

        self._test_contributions(
            shaker,
            [('FFFF', 1), ('XDDD', 1)],
            [derivatives_e.convert_frequency_from_string('1500nm')],
            ['FF_FF__0_0', 'F_FFF__0_0', 'F_F_FF__1_0', 'F_F_FF__0_1'],
            [
                'FF_FF__1_1', 'FF_FF__2_0', 'FF_FF__0_2',
                'F_FFF__1_1', 'F_FFF__2_0', 'F_FFF__0_2',
                'F_F_F_F__1_1', 'F_F_F_F__2_0', 'F_F_F_F__0_2'
            ],
            is_zpva=True)

        self._test_contributions(
            shaker,
            [('XDDD', 2)],
            [derivatives_e.convert_frequency_from_string('1500nm')],
            [
                'FF_FF__0_0', 'FF_FF__1_1', 'FF_FF__2_0', 'FF_FF__0_2',
                'F_FFF__0_0', 'F_FFF__1_1', 'F_FFF__2_0', 'F_FFF__0_2',
                'F_F_FF__1_0', 'F_F_FF__0_1',
                'F_F_F_F__1_1', 'F_F_F_F__2_0', 'F_F_F_F__0_2'
            ], is_zpva=True)

    def test_shaking_save_and_load(self):
        """Test that we are able to save and load vibrational contributions"""

        df = chemistry_datafile.ChemistryDataFile()

        with open(self.datafile) as f:
            df.read(f)

        shaker = shaking.Shaker(datafile=df)

        must_be_in = ['F_F__0_0', 'F_F__1_1', 'F_F__2_0', 'F_F__0_2']
        vibs = self._test_contributions(
            shaker,
            [('FF', 2), ('dD', 2)],
            [0.02, 0.04],
            must_be_in)

        # test saving
        shaking.save_vibrational_contributions(self.datafile, vibs)

        # test loading:
        nvibs = shaking.load_vibrational_contributions(self.datafile, df.spacial_dof)

        self.assertIn('FF', nvibs)
        self.assertIn('dD', nvibs)

        for i in ['FF', 'dD']:
            c = vibs[i]
            cx = nvibs[i]
            for j in c.vibrational_contributions:
                self.assertIn(j, cx.vibrational_contributions)
                for freq in c.vibrational_contributions[j]:
                    self.assertIn(freq, cx.vibrational_contributions[j])
                    self.assertTensorsAlmostEqual(
                        c.vibrational_contributions[j][freq], cx.vibrational_contributions[j][freq])

                for freq in c.total_zpva:
                    self.assertIn(freq, cx.total_zpva)
                    self.assertTensorsAlmostEqual(
                        c.total_zpva[freq], cx.total_zpva[freq])
                    self.assertIn(freq, cx.total_pv)
                    self.assertTensorsAlmostEqual(
                        c.total_pv[freq], cx.total_pv[freq])
                    self.assertIn(freq, cx.total_vibrational)
                    self.assertTensorsAlmostEqual(
                        c.total_vibrational[freq], cx.total_vibrational[freq])

    def test_nachos_shake(self):
        """Test the shake command"""

        must_be_in = ['F_F__0_0', 'F_F__1_1', 'F_F__2_0', 'F_F__0_2']
        d = ['FF', 'dD']
        notd = ['F', 'FFF', 'dDF', 'XDD']
        freqs = [0.02, 0.04]

        # process
        process = self.run_python_script(
            'nachos/shake.py', ['-d', self.datafile, '-O', ';'.join(d), '-f', ';'.join(str(a) for a in freqs)],
            out_pipe=subprocess.PIPE,
            err_pipe=subprocess.PIPE)

        stdout_t, stderr_t = process.communicate()

        self.assertEqual(len(stderr_t), 0, msg=stderr_t.decode())
        self.assertEqual(len(stdout_t), 0, msg=stdout_t.decode())

        nvibs = shaking.load_vibrational_contributions(self.datafile, 15)

        for i in d:
            self.assertIn(i, nvibs)
            cx = nvibs[i]
            for j in must_be_in:
                self.assertIn(j, cx.vibrational_contributions)

                if 'D' in i:
                    for f in freqs:
                        self.assertIn(f, cx.vibrational_contributions[j])
                else:
                    self.assertIn('static', cx.vibrational_contributions[j])
        for i in notd:
            self.assertNotIn(i, nvibs)
