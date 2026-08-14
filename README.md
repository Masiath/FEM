# FEM

Finite Element Method (FEM) analysis projects using [NGSolve](https://ngsolve.org/).

## Overview

This repository is a collection of NGSolve-based FEM simulation projects. Each project lives in its own subfolder under `projects/` and is self-contained (own code, data, results, and figures).

## Repository structure

```
FEM/
├── projects/
│   ├── CystatinC/
│   │   ├── src/          # NGSolve simulation scripts (mesh generation, solvers, sweeps)
│   │   ├── data/         # Input parameters, material data, refractive index tables
│   │   ├── results/      # Output data (CSV/JSON) from simulation runs
│   │   ├── figures/      # Plots generated from results
│   │   ├── notebooks/    # Exploratory/analysis notebooks
│   │   └── README.md     # Project-specific overview, methodology, and results
│   └── .../              # Additional projects follow the same layout
├── requirements.txt
└── README.md
```

Each project's own README documents its purpose, methodology, and results — this top-level README stays project-agnostic.

## Setup

```bash
python -m venv venv
source venv/bin/activate        # or venv\Scripts\activate on Windows
pip install -r requirements.txt
```

Requires NGSolve (`pip install ngsolve`) and standard scientific Python packages (numpy, scipy, matplotlib, plotly).

## Projects

- **CystatinC** — FEM simulation project (see `projects/CystatinC/README.md` for details).

## Adding a new project

1. Create a new folder under `projects/` with a descriptive name.
2. Follow the `src/ data/ results/ figures/ notebooks/` layout above.
3. Add a project-level `README.md` describing what it does.

## Figure style

For consistency across projects: individual plots only (no combined subplots), axis label font size 18, axis tick/value font size 22, title font size 22.

## License

[MIT](LICENSE)
