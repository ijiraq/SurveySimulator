"""Smoke tests for RosePlot (Agg backend; no display)."""
from __future__ import annotations

import pathlib
import tempfile
import unittest

import matplotlib

matplotlib.use('Agg')

from astropy import units as u
from astropy.table import QTable
from astropy.time import Time

from ossssim.plotter import RosePlot, _as_time, _wedge_width_from_area

script_directory = pathlib.Path(__file__).parent.resolve()
CFEPS = script_directory / 'data' / 'Surveys' / 'CFEPS'


class TestRosePlotHelpers(unittest.TestCase):
    def test_as_time_accepts_time_and_jd(self):
        t = Time('2020-01-01')
        self.assertIsInstance(_as_time(t), Time)
        self.assertIsInstance(_as_time(t.jd), Time)

    def test_as_time_rejects_angle_quantity(self):
        with self.assertRaises(TypeError):
            _as_time(30 * u.deg)

    def test_wedge_width_from_area(self):
        w = _wedge_width_from_area(4.0)
        self.assertAlmostEqual(w.to(u.deg).value, 2.0)


class TestRosePlotSmoke(unittest.TestCase):
    def test_construct_and_savefig(self):
        plot = RosePlot(Time('2020-01-01'))
        plot.add_scale_rings([20, 40])
        with tempfile.NamedTemporaryFile(suffix='.png') as tmp:
            plot.savefig(tmp.name)
            self.assertGreater(pathlib.Path(tmp.name).stat().st_size, 0)

    def test_add_model_cartesian(self):
        table = QTable({
            'x': [30.0, 40.0] * u.au,
            'y': [0.0, 10.0] * u.au,
            'z': [0.0, 0.0] * u.au,
        })

        class _Tiny:
            def __init__(self, table):
                self.table = table

        plot = RosePlot(Time('2020-01-01'))
        plot.add_model(_Tiny(table), ms=2)
        with tempfile.NamedTemporaryFile(suffix='.png') as tmp:
            plot.savefig(tmp.name)

    def test_add_pointings_cfeps(self):
        if not CFEPS.is_dir():
            self.skipTest(f'CFEPS fixture missing: {CFEPS}')
        try:
            import ossssimlib  # noqa: F401
        except ImportError:
            self.skipTest('ossssimlib not importable')
        plot = RosePlot(Time('2020-01-01'))
        n = plot.add_pointings(CFEPS)
        self.assertGreater(n, 0)


if __name__ == '__main__':
    unittest.main()
