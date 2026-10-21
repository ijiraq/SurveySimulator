"""Geocentric (code 500) sky positions near an equinox.

Earth's ICRF z is below 0.1 r within ~15 deg of the equinoxes. RADECeclXV
once took such observer vectors for ecliptic and rotated them again, which
moved a 45 au object by ~150" and turned its direction of motion by ~27 deg.
"""
from __future__ import annotations

import pathlib
import tempfile
import unittest

from astropy import units as u

from ossssim import OSSSSim
from ossssim.color import PhotSpec

# 2019 QF111 (DEEP VI Table 2), barycentric ecliptic elements at JD 2458725.5,
# seen from the geocentre on JD 2458754.59448 (2019-09-28, a week after the
# September equinox). Expected RA/Dec: the object rotated to ICRF minus the
# geocentre in ICRF (JPL DE ephemeris through astropy), no light time.
ELEMENTS = dict(a=44.3076856 * u.au, e=0.0792533, inc=3.3757953 * u.deg,
                node=33.8367 * u.deg, peri=186.0071366 * u.deg, M=124.47668073 * u.deg,
                epoch=2458725.5 * u.day)
JD = 2458754.59448
RA, DEC = 352.870225, -5.607840
HALF = 30.0 / 3600.0     # half-width of the field [deg]

EFF = """\
rate_cut= 0.0 20.0 0.0 180.0
track_frac= 1.0 99.0 0.0
filter= r
rates= 0.00 20.00
function= lookup
lookup_param= 0.0 1.0
lookup_param= 40.0 1.0
mag_lim= 40.0
"""


class GeocentricEquinoxTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = pathlib.Path(self.tmp.name) / 'Equinox'
        root.mkdir()
        corners = [(-HALF, -HALF), (HALF, -HALF), (HALF, HALF), (-HALF, HALF)]
        lines = [f'poly 4 {RA:.6f} {DEC:.6f} {JD:.5f} 1.0 500 field.eff']
        lines += [f'{x:+.6f} {y:+.6f}' for x, y in corners]
        (root / 'pointings.list').write_text('\n'.join(lines) + '\n')
        (root / 'field.eff').write_text(EFF)
        self.root = root

    def tearDown(self):
        self.tmp.cleanup()

    def test_object_lands_in_a_one_arcminute_field(self):
        sim = OSSSSim(str(self.root), seed=11)
        row = dict(ELEMENTS, H=6.0 * u.mag, comp='default')
        out = sim.simulate(row, colors=PhotSpec(), model_band='r')
        self.assertGreaterEqual(int(out['flag']), 4)
        self.assertAlmostEqual(out['RA'].to(u.deg).value, RA, delta=1.0 / 3600)
        self.assertAlmostEqual(out['DEC'].to(u.deg).value, DEC, delta=1.0 / 3600)


if __name__ == '__main__':
    unittest.main()
