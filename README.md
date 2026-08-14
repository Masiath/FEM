# FEM

Finite Element Method (FEM) analysis using [NGSolve](https://ngsolve.org/) for validating surface plasmon resonance (SPR) biosensor designs against transfer matrix method (TMM) results.

## Overview

This repository contains NGSolve-based FEM simulations used to independently verify the optical response of multilayer SPR biosensor stacks (e.g. prism/metal/dielectric/analyte configurations) originally modeled via TMM. The FEM approach solves the full-wave electromagnetic field problem on a 2D/axisymmetric mesh, providing a cross-check on resonance angle and sensitivity predictions.

## Status

Actively in progress. Current validated results:

| Quantity | Value |
|---|---|
| SPR resonance angle (θ_SPR) | 84.7375° |
| Sensitivity | ≈ 547–548 °/RIU |

## Project structure

```
FEM/
├── src/            # NGSolve simulation scripts (mesh generation, solvers, sweeps)
├── data/           # Input parameters, material dispersion data, refractive index tables
├── results/        # Output data (CSV/JSON) from simulation runs
├── figures/        # Publication-quality plots generated from results
├── notebooks/      # Exploratory/analysis notebooks
├── requirements.txt
└── README.md
```

## Setup

```bash
python -m venv venv
source venv/bin/activate        # or venv\Scripts\activate on Windows
pip install -r requirements.txt
```

Requires NGSolve (`pip install ngsolve`) and standard scientific Python packages (numpy, scipy, matplotlib, plotly).

## Usage

Simulation scripts live in `src/`. Each script generates the mesh for a given layer stack, solves the Helmholtz/Maxwell system at a sweep of incidence angles, and extracts the reflectance curve to locate θ_SPR. Results are written to `results/`, and figures are generated from `figures/` scripts using Plotly (per-figure files, no combined subplots).

## Notes

- Figures follow a fixed style: individual plots (no subplots), axis label font size 18, axis tick/value font size 22, title font size 22.
- FEM results are cross-validated against a TMM/brute-force optimization pipeline maintained separately.

## License

TBD
