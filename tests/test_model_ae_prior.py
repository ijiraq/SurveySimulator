"""Unit tests for JWST model-conditioned p(a,e|r,i) debiasing helpers."""
from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "JWST" / "scripts" / "grid_bias.py"
L7 = ROOT / "F95" / "tests" / "Models" / "L7model-3.0-9.0"

spec = importlib.util.spec_from_file_location("grid_bias_model_ae", HELPER)
grid_bias = importlib.util.module_from_spec(spec)
sys.modules["grid_bias_model_ae"] = grid_bias
spec.loader.exec_module(grid_bias)


def _write_mini_model(path: Path) -> None:
    """Tiny L7-style catalog with cold-like and hot-like rows."""
    lines = [
        "# Epoch of elements: JD = 2453157.50000",
        "# Longitude of Neptune: lambdaN = 5.489",
        "#",
        "#   a      e     i     node    peri    M       H      dist     comp      j  k",
    ]
    # cold classical-ish at r~44, i~2
    for k in range(40):
        a = 44.0 + 0.01 * k
        e = 0.04
        i = 1.5 + 0.02 * k
        dist = 43.5 + 0.02 * k
        lines.append(
            f"  {a:.3f}  {e:.3f}  {i:.3f}  10.0  20.0  30.0  8.0  {dist:.2f}   c  0  0"
        )
    # hot classical-ish at r~44, i~20
    for k in range(40):
        a = 42.0 + 0.05 * k
        e = 0.15
        i = 18.0 + 0.1 * k
        dist = 42.0 + 0.05 * k
        lines.append(
            f"  {a:.3f}  {e:.3f}  {i:.3f}  10.0  20.0  30.0  8.0  {dist:.2f}   h  0  0"
        )
    # resonant at large r / moderate i (for expand test)
    for k in range(10):
        a = 55.0
        e = 0.25
        i = 10.0
        dist = 60.0 + 0.1 * k
        lines.append(
            f"  {a:.3f}  {e:.3f}  {i:.3f}  10.0  20.0  30.0  8.0  {dist:.2f}   resonant  5  1"
        )
    path.write_text("\n".join(lines) + "\n")


class OrbitModelCatalogTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.model_path = Path(self.tmp.name) / "mini_l7.txt"
        _write_mini_model(self.model_path)
        self.model = grid_bias.OrbitModelCatalog.from_l7(self.model_path)

    def tearDown(self):
        self.tmp.cleanup()

    def test_load_counts(self):
        self.assertEqual(len(self.model), 90)

    def test_select_mixes_components_by_i(self):
        cold = self.model.select(43.0, 45.0, 1.0, 3.0)
        hot = self.model.select(41.0, 45.0, 17.0, 25.0)
        self.assertGreater(len(cold), 10)
        self.assertGreater(len(hot), 10)
        self.assertTrue(all(c == "c" for c in cold.comp))
        self.assertTrue(all(c == "h" for c in hot.comp))

    def test_low_i_window_is_not_forced_to_one_component_label(self):
        # Real L7 mixes c and f at low i; here only c exists at low i.
        # Document that we do not filter by Sample A cold/hot tags.
        self.assertEqual(grid_bias.JWST_SAMPLE_A.bias_method, "model_ae")
        self.assertEqual(grid_bias.N26_HELIOSTACK.bias_method, "aq_grid")

    def test_sample_ae_respects_r(self):
        rng = np.random.default_rng(0)
        cold = self.model.select(43.0, 45.0, 1.0, 3.5)
        for _ in range(50):
            a, e, comp = cold.sample_ae(rng, r_au=44.0)
            q = a * (1.0 - e)
            q_ap = a * (1.0 + e)
            self.assertLessEqual(q, 44.0 + 1e-8)
            self.assertGreaterEqual(q_ap, 44.0 - 1e-8)
            self.assertEqual(comp, "c")

    def test_select_expanding_grows_sparse_window(self):
        # Tiny window around resonant objects; min_n forces expand.
        sub, dr, di = self.model.select_expanding(
            60.0, 60.2, 10.0, 10.2, min_n=8, max_expand=6
        )
        self.assertGreaterEqual(len(sub), 8)
        self.assertGreater(dr, 0.1)

    def test_rih_cell_key(self):
        key = grid_bias.rih_cell_key(46.20, 1.85, 8.12)
        self.assertEqual(key[0], 46.0)
        self.assertEqual(key[1], 1.0)
        bounds = grid_bias.rih_bounds_from_key(key)
        self.assertAlmostEqual(bounds["r"][1] - bounds["r"][0], grid_bias.R_STEP)
        self.assertAlmostEqual(bounds["i"][1] - bounds["i"][0], grid_bias.I_STEP)

    def test_jwst_detections_use_rih_cells(self):
        dets = grid_bias.load_detections(
            ROOT / "JWST" / "data" / "jwst_sampleA.csv",
            grid_bias.JWST_SAMPLE_A,
        )
        self.assertEqual(len(dets), 20)
        # JPB04: r=46.20, i=1.85 → cell (46.0, 1.0, ...)
        jpb04 = next(d for d in dets if d["name"] == "JPB04")
        self.assertEqual(jpb04["cell"][0], 46.0)
        self.assertEqual(jpb04["cell"][1], 1.0)
        self.assertEqual(len(jpb04["cell"]), 3)

    def test_n26_detections_still_use_aq_cells(self):
        n26 = grid_bias.N26_HELIOSTACK
        path = ROOT / "N26" / "data" / "n26_detections.csv"
        if not path.exists():
            self.skipTest("N26 catalog missing")
        dets = grid_bias.load_detections(path, n26)
        self.assertEqual(len(dets[0]["cell"]), 4)

    def test_sample_aimed_elements_at_i_pins_r_and_i(self):
        rng = np.random.default_rng(1)
        jpl = ROOT / "JWST" / "characterization" / "epoch1" / "JWST.csv"
        obs = grid_bias.parse_jpl_horizons_icrf(
            jpl, grid_bias.JWST_SAMPLE_A.epoch_jd[0]
        )
        a, e, i_deg, r_au = 44.0, 0.05, 2.5, 44.0
        el = grid_bias.sample_aimed_elements_at_i(
            a, e, i_deg, obs, rng, r_au, survey=grid_bias.JWST_SAMPLE_A
        )
        self.assertIsNotNone(el)
        inc, node, peri, M = el
        self.assertAlmostEqual(inc, i_deg, places=6)
        # Recovered heliocentric distance should match the pin.
        xyz = grid_bias.ecliptic_xyz_from_elements(a, e, inc, node, peri, M)
        r_got = float(np.linalg.norm(xyz))
        self.assertAlmostEqual(r_got, r_au, places=4)


@unittest.skipUnless(L7.exists(), "L7 model fixture not present")
class L7ModelPriorIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = grid_bias.OrbitModelCatalog.from_l7(L7)

    def test_l7_low_i_mixes_classical_components(self):
        sub, _dr, _di = self.model.select_expanding(45.0, 46.0, 1.0, 2.0, min_n=50)
        fracs = sub.component_fractions()
        # At low i near 45 AU, L7 is dominated by classical (c/f), not a
        # single forced Sample-A cold/hot label.
        classical = fracs.get("c", 0.0) + fracs.get("f", 0.0)
        self.assertGreater(classical, 0.8)
        self.assertGreater(len(fracs), 1)

    def test_l7_high_i_is_hot_or_resonant(self):
        sub, _dr, _di = self.model.select_expanding(41.0, 42.0, 26.0, 27.0, min_n=50)
        fracs = sub.component_fractions()
        self.assertGreater(fracs.get("h", 0.0) + fracs.get("resonant", 0.0), 0.9)
        self.assertLess(fracs.get("c", 0.0) + fracs.get("f", 0.0), 0.05)


if __name__ == "__main__":
    unittest.main()
