"""Tests for SurveyCharacterization area / fill / depth coverage API."""
import pathlib
import unittest

import numpy as np

from ossssim import SurveyCharacterization

script_directory = pathlib.Path(__file__).parent.resolve()
CFEPS = script_directory / 'data' / 'Surveys' / 'CFEPS'
COSDEC = script_directory.parent / 'F95' / 'tests' / 'data' / 'cosdec_survey'
COSDEC_GC = script_directory.parent / 'F95' / 'tests' / 'data' / 'cosdec_gc'


class SurveyCharacterizationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not CFEPS.is_dir():
            raise unittest.SkipTest(f'CFEPS fixture missing: {CFEPS}')
        try:
            cls.survey = SurveyCharacterization.from_directory(str(CFEPS))
        except ImportError:
            raise unittest.SkipTest('ossssimlib not importable')

    def test_loads_pointings(self):
        self.assertGreater(len(self.survey), 0)
        self.assertIn('L3f-smooth', self.survey.keys())
        self.assertIn('L3h-smooth', self.survey.pointings)

    def test_duplicate_eff_keys_disambiguated(self):
        # L3w has two pointings sharing L3w-smooth.eff
        keys = self.survey.keys()
        self.assertTrue(
            any(k.startswith('L3w-smooth') for k in keys),
            msg=f'expected L3w-smooth keys among {keys[:10]}...',
        )
        l3w = [k for k in keys if k.startswith('L3w-smooth')]
        self.assertGreaterEqual(len(l3w), 2)
        self.assertIn('L3w-smooth#0', l3w)
        self.assertIn('L3w-smooth#1', l3w)

    def test_area_and_fill(self):
        p = self.survey['L3f-smooth']
        self.assertGreater(p.area_deg2, 0.0)
        self.assertGreater(p.fill_factor, 0.0)
        self.assertLessEqual(p.fill_factor, 1.0)
        self.assertAlmostEqual(p.fill_factor, 0.80, places=2)

    def test_efficiency_and_effective_area(self):
        p = self.survey['L3h-smooth']
        eta = p.efficiency(23.0)
        self.assertGreater(eta, 0.0)
        self.assertLessEqual(eta, 1.0)
        a_eff = p.effective_area(23.0)
        self.assertAlmostEqual(
            a_eff, p.area_deg2 * p.fill_factor * eta, places=10
        )

    def test_coverage_vs_magnitude_decreases(self):
        mags = np.arange(21.0, 26.01, 0.5)
        m, area = self.survey.coverage_vs_magnitude(mags)
        self.assertEqual(list(m), list(mags))
        self.assertTrue(np.all(area >= 0.0))
        # Fainter end should not be larger than bright end for survey coverage
        self.assertLessEqual(area[-1], area[0] + 1e-9)
        # Check consistency with single-pointing sum at one mag
        mid = 23.0
        _, a_mid = self.survey.coverage_vs_magnitude([mid])
        expected = sum(p.effective_area(mid) for p in self.survey)
        self.assertAlmostEqual(float(a_mid[0]), float(expected), places=8)

    def test_subset_pointing_ids(self):
        mags = np.array([22.0, 24.0])
        _, a_all = self.survey.coverage_vs_magnitude(mags)
        _, a_one = self.survey.coverage_vs_magnitude(
            mags, pointing_ids=['L3f-smooth']
        )
        self.assertTrue(np.all(a_one <= a_all + 1e-9))
        p = self.survey['L3f-smooth']
        np.testing.assert_allclose(
            a_one, p.effective_area(mags), rtol=1e-10
        )


class CosDecFootprintTest(unittest.TestCase):
    """1 deg x 1 deg footprints must have the same area at every Dec."""

    @classmethod
    def setUpClass(cls):
        if not COSDEC.is_dir():
            raise unittest.SkipTest(f'cosdec fixture missing: {COSDEC}')
        try:
            cls.survey = SurveyCharacterization.from_directory(str(COSDEC))
        except ImportError:
            raise unittest.SkipTest('ossssimlib not importable')

    def test_square_area_independent_of_dec(self):
        # Index 6 is the MegaPrime ears footprint, 8 the 8x1 deg strip and
        # 9 the slanted OSSOS-like block.
        squares = [p for i, p in enumerate(self.survey) if i not in (6, 8, 9)]
        self.assertEqual(len(squares), 7)
        areas = np.array([p.area_deg2 for p in squares])
        decs = np.degrees([p.dec for p in squares])
        self.assertLess(abs(decs).max(), 90.0)
        self.assertGreater(np.ptp(decs), 100.0)
        np.testing.assert_allclose(areas, 1.0, rtol=1e-3)

    def test_strip_area(self):
        self.assertAlmostEqual(self.survey[8].area_deg2 / 8.0, 1.0, delta=1e-3)


class GreatCircleEdgesTest(unittest.TestCase):
    """survey.conf 'edges = great_circle' keeps long edges undensified."""

    def test_strip_area_great_circle(self):
        if not COSDEC_GC.is_dir():
            raise unittest.SkipTest(f'fixture missing: {COSDEC_GC}')
        try:
            survey = SurveyCharacterization.from_directory(str(COSDEC_GC))
        except ImportError:
            raise unittest.SkipTest('ossssimlib not importable')
        self.assertEqual(len(survey), 1)
        # Great-circle top and bottom edges bulge poleward: ~0.65% more area
        # than the 8 deg^2 constant-Dec strip.
        self.assertGreater(survey[0].area_deg2, 8.03)
        self.assertLess(survey[0].area_deg2, 8.08)


if __name__ == '__main__':
    unittest.main()
