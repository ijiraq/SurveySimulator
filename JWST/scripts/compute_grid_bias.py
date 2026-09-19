#!/usr/bin/env python3
"""Grid-cell debiasing for JWST Sample A (Kavelaars et al. 2022 ac2c72 method)."""
from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path

import numpy as np
from astropy import units as u
from ossssim import OSSSSim
from ossssim.color import PhotSpec

A_STEP = 0.2
Q_STEP = 0.2
SI_STEP = 0.001
H_STEP = 0.1
TARGET_DETECTIONS = 5000
EPOCH_JD = [2459969.5, 2459974.5, 2459978.5]
FIELD_RA = 209.3875
FIELD_DEC = -10.865278
H_COLOR_OFFSET = 1.0


def laplace_inclination(a_au: float) -> float:
    return 1.759 + 0.0321 * (a_au - 41.8)


def laplace_node(a_au: float) -> float:
    return 90.0 - 0.5 * (a_au - 43.0)


def compute_ifree(i_deg: float, omega_deg: float, a_au: float) -> float:
    ip = laplace_inclination(a_au)
    om_lp = laplace_node(a_au)
    cos_ifree = (
        math.cos(math.radians(i_deg)) * math.cos(math.radians(ip))
        + math.sin(math.radians(i_deg)) * math.sin(math.radians(ip))
        * math.cos(math.radians(omega_deg - om_lp))
    )
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_ifree))))


def apparent_to_Hr(m_f150w2: float, d_au: float, phase: float = 0.35) -> float:
    return m_f150w2 + H_COLOR_OFFSET - 10.0 * math.log10(d_au) + phase


def cell_index(value: float, step: float) -> float:
    return math.floor(value / step) * step


def cell_key(a: float, q: float, sin_ifree: float, hx: float) -> tuple:
    return (
        round(cell_index(a, A_STEP), 6),
        round(cell_index(q, Q_STEP), 6),
        round(cell_index(sin_ifree, SI_STEP), 6),
        round(cell_index(hx, H_STEP), 6),
    )


def load_detections(path: Path) -> list[dict]:
    rows = []
    with path.open() as fh:
        for row in csv.DictReader(fh):
            a, e, i = float(row["a"]), float(row["e"]), float(row["i"])
            d = float(row["d_bary"])
            hx = apparent_to_Hr(float(row["m_f150w2"]), d)
            q = a * (1.0 - e)
            ifree = compute_ifree(i, 0.0, a)
            rows.append({**row, "a": a, "e": e, "i": i, "d_bary": d, "q": q,
                         "Hx": hx, "ifree": ifree,
                         "sin_ifree": math.sin(math.radians(ifree)),
                         "cell": cell_key(a, q, math.sin(math.radians(ifree)), hx)})
    return rows


def load_bias_cache(path: Path) -> dict:
    if not path.exists():
        return {}
    out = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            key = tuple(float(row[k]) for k in ("a_bin", "q_bin", "si_bin", "h_bin"))
            out[key] = (float(row["bias"]), int(row["n_drawn"]))
    return out


def save_bias_cache(path: Path, cache: dict) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["a_bin", "q_bin", "si_bin", "h_bin", "bias", "n_drawn"])
        for key, (bias, n_drawn) in sorted(cache.items()):
            w.writerow([*key, bias, n_drawn])


def setup_pointings(char_root: Path) -> None:
    # Implant search area 1.6x1.6 deg; active mosaic 0.05 deg2
    ff = 0.05 / (1.6 * 1.6)
    for idx, jd in enumerate(EPOCH_JD, start=1):
        text = (
            f"# JWST Sample A epoch {idx}\n"
            f"1.6 1.6 {FIELD_RA} {FIELD_DEC} {jd} {ff:.5f} JWST.csv JWST_sampleA.eff\n"
        )
        (char_root / f"epoch{idx}" / "pointings.list").write_text(text)


class JWSTSimulator:
    def __init__(self, char_root: Path, seed: int = 42):
        os.chdir("/arc/home/jkavelaars/Develop/SurveySimulator")
        setup_pointings(char_root)
        self.sims = [OSSSSim(str(char_root / f"epoch{i}"), seed=seed + i) for i in (1, 2, 3)]
        self.colors = PhotSpec()

    def detected_sample_a(self, a, e, inc, node, peri, M, H) -> bool:
        base = dict(a=a * u.au, e=e, inc=inc * u.deg, node=node * u.deg, peri=peri * u.deg,
                    M=M * u.deg, H=H * u.mag, comp="default")
        for sim, jd in zip(self.sims, EPOCH_JD):
            r = sim.simulate({**base, "epoch": jd * u.day}, colors=self.colors, model_band="r")
            if r["flag"] < 4:
                return False
        return True


def compute_cell_bias(sim: JWSTSimulator, cell_bounds: dict, seed: int, target: int) -> tuple[float, int]:
    rng = np.random.default_rng(seed)
    a0, a1 = cell_bounds["a"]
    q0, q1 = cell_bounds["q"]
    si0, si1 = cell_bounds["sin_ifree"]
    h0, h1 = cell_bounds["Hx"]

    n_detected = 0
    n_drawn = 0
    max_draws = max(target * 200000, 500000)
    while n_detected < target and n_drawn < max_draws:
        a = rng.uniform(a0, a1)
        q = rng.uniform(q0, min(q1, a * 0.98))
        if q <= 0:
            continue
        e = 1.0 - q / a
        sin_ifree = rng.uniform(si0, si1)
        ifree = math.degrees(math.asin(max(-1.0, min(1.0, sin_ifree))))
        H = rng.uniform(h0, h1)
        node, peri, M = rng.uniform(0, 360, size=3)
        n_drawn += 1
        if sim.detected_sample_a(a, e, max(ifree, 0.05), node, peri, M, H):
            n_detected += 1
        if n_drawn % 50000 == 0:
            print(f"    ... {n_drawn} draws, {n_detected}/{target} detections", flush=True)
    if n_detected < target:
        raise RuntimeError(f"Only {n_detected}/{target} after {n_drawn} draws")
    return n_detected / n_drawn, n_drawn


def bounds_from_key(key: tuple) -> dict:
    a0, q0, si0, h0 = key
    return {"a": (a0, a0 + A_STEP), "q": (q0, q0 + Q_STEP),
            "sin_ifree": (max(0.0, si0), si0 + SI_STEP), "Hx": (h0, h0 + H_STEP)}


def write_detections_full(out_path: Path, detections: list[dict]) -> None:
    header = f"""# File: JWST-free-cla_m.detections-full
#
# Grid debiasing ac2c72; Eduardo et al. 2026 Sample A (20 objects)
# H_r from m_F150W2 + 1.0 - 10log10(d) + 0.35
#
# Grid size:
# h_step:  {H_STEP}
# a_step:  {A_STEP}
# q_step:  {Q_STEP}
# si_step: {SI_STEP}
#
"""
    cols = ("cl p j k sh object mag e_mag Filt Hsur dist e_dist Nobs time av_xres av_yres "
            "max_x max_y a e_a e i e_i Omega e_Omega omega e_omega tperi e_tperi "
            "RAdeg DEdeg JD rate MPC ifree Omfree omfree Hx comp bias")
    lines = [header, cols]
    for d in detections:
        lines.append(
            f"cla m -1 -1 S {d['name']:7s} {d['Hx']:.2f} 0.100 r {d['Hx']:.2f} {d['d_bary']:.3f} 0.100 "
            f"3 0.0000 0.083 0.073 0.311 0.343 {d['a']:11.6f} 0.1012 {d['e']:.6f} 0.001009 "
            f"{d['i']:6.3f} 0.100 0.000 0.100 0.000 0.100 0.000 0.100 0.000 0.100 "
            f"{FIELD_RA:.3f} {FIELD_DEC:.3f} {EPOCH_JD[1]:.5f} 0.40 {d['name']:7s} {d['ifree']:6.3f} 0.000 0.000 "
            f"{d['Hx']:.2f} {d['comp']} {d['bias']:.7f}"
        )
    out_path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--target", type=int, default=TARGET_DETECTIONS)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    root = Path(args.root)
    detections = load_detections(root / "data/jwst_sampleA.csv")
    cache_path = root / "bias_grid.csv"
    cache = load_bias_cache(cache_path)
    sim = JWSTSimulator(root / "characterization", seed=args.seed)
    cells = sorted({d["cell"] for d in detections})
    print(f"{len(cells)} cells, target={args.target}/cell")

    for idx, key in enumerate(cells):
        if key in cache:
            print(f"cell {idx+1}/{len(cells)} {key}: cached {cache[key][0]:.4g}")
            continue
        print(f"cell {idx+1}/{len(cells)} {key}:")
        bias, n_drawn = compute_cell_bias(sim, bounds_from_key(key), args.seed + idx, args.target)
        cache[key] = (bias, n_drawn)
        print(f"  bias={bias:.6g} n_drawn={n_drawn}")
        save_bias_cache(cache_path, cache)

    for d in detections:
        d["bias"] = cache[d["cell"]][0]
    out = root / "JWST-free-cla_m.detections-full"
    write_detections_full(out, detections)
    print(f"Wrote {out}; sum 1/bias = {sum(1/d['bias'] for d in detections):.1f}")


if __name__ == "__main__":
    main()
