# OSSOS Models 1.0 — ModelUsed tables

Orbit tables used as the empirical \(p(a,e\mid r,i)\) prior for JWST
Sample A grid debiasing (`bias_method=model_ae`).

## Provenance

Copied from the public [OSSOS/OSSOS_Models](https://github.com/OSSOS/OSSOS_Models)
component branches (`ModelUsed-check-8.66.dat` on each of Classical,
Detached, Inner, Plutinos, Scattering, Twotinos). These are the OSSOS 1.0
nominal populations calibrated to the OSSOS++ sample
(Bannister et al. 2018, ApJS 236, 18).

| File | Source branch / path |
|------|----------------------|
| `Classical-ModelUsed.dat` | `Classical/Classical/ModelUsed-check-8.66.dat` |
| `Detached-ModelUsed.dat` | `Detached/Detached/ModelUsed-check-8.66.dat` |
| `Inner-ModelUsed.dat` | `Inner/Inner/ModelUsed-check-8.66.dat` |
| `Plutinos-ModelUsed.dat` | `Plutinos/Plutinos/ModelUsed-check-8.66.dat` |
| `Scattering-ModelUsed.dat` | `Scattering/Scattering/ModelUsed-check-8.66.dat` |
| `Twotinos-ModelUsed.dat` | `Twotinos/Twotinos/ModelUsed-check-8.66.dat` |

## Format

OSSOS `ModelUsed.dat` columns:

```text
a  e  i  Omega  omega  M  H  epoch  dist  comment
```

`ossssim.grid_bias.OrbitModelCatalog.from_path` loads this directory and
tags each row with the filename population (`Classical`, `Scattering`, …).

## Configuration

```bash
# default: this directory (all components)
python JWST/scripts/compute_grid_bias.py

# single component or an alternate table
python JWST/scripts/compute_grid_bias.py --model Models/OSSOS/Classical-ModelUsed.dat

# legacy CFEPS L7 fixture still works
python JWST/scripts/compute_grid_bias.py --model F95/tests/Models/L7model-3.0-9.0
```
