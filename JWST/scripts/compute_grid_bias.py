#!/usr/bin/env python3
"""Grid-cell debiasing for JWST Sample A (model-conditioned a/e prior)."""
from __future__ import annotations

import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from ossssim.grid_bias import H_STEP, I_STEP, JWST_SAMPLE_A, R_STEP
from ossssim.grid_bias_run import main

_HEADER = f"""# File: JWST-free-cla_m.detections-full
#
# Grid debiasing ac2c72; Eduardo et al. 2026 Sample A (20 objects)
# Bias method: model_ae — cells in (r, ecliptic i, H)
# (a, e) drawn from OSSOS Models 1.0 ModelUsed tables (Models/OSSOS/;
# override with --model FILE|DIR) conditioned on the cell's (r, i) window.
# All dynamical components are kept (mixture = model prior).
# Catalog cold/hot tags are not used to filter the prior.
# H_r from m_F150W2 + 1.0 - 5log10(r Δ) + 2.5log10(Bowell Φ), r=Δ=d_bary, G=-0.12
# Extra after MPC: ifree Omfree omfree (Laplace-free; Omfree=omfree=0), Hx, comp, bias
#
# Grid size:
# h_step:  {H_STEP}
# r_step:  {R_STEP}
# i_step:  {I_STEP}
#
"""


if __name__ == "__main__":
    main(
        JWST_SAMPLE_A,
        default_root=Path(__file__).resolve().parents[1],
        extra_header=_HEADER,
    )
