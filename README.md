# Symmetry-Based Magnon Band Analysis Toolkit

This repository contains an in-house computational workflow for symmetry-based analysis of magnetic materials, with a focus on magnon band structures and symmetry-protected degeneracies.

Starting from a crystal structure (`POSCAR_JULIA`) and magnetic configuration (`INFO_JULIA` / `INFO_FORT`), the pipeline can:

- Identify the full symmetry group of the system, including space groups, magnetic space groups, and spin groups.
- Construct little-group representations at arbitrary `k` points.
- Decompose magnon bands into irreducible representations.
- Analyze symmetry-protected band degeneracies.
- Generate symmetry-constrained neighbor information for subsequent magnon calculations.

The generated neighbor data are then passed to the magnon solver, which computes the corresponding band structure. In this way, the toolkit connects symmetry analysis directly to physical observables:

`Structure -> Symmetry -> Hamiltonian -> Band -> Representation -> Degeneracy`

## Why This Project

This project is intended as a research-oriented software portfolio piece. It demonstrates that the underlying algorithms, data processing workflow, and numerical pipeline were implemented in code rather than described only at a conceptual level.

## Repository Structure

- `symmetry-analyzer/`: Julia code for symmetry detection, representation analysis, and supporting databases.
- `magnon-solver/`: Fortran-based magnon solver and Julia interface code.
- `demo/`: runnable examples and outputs prepared for quick testing by collaborators, interviewers, or reviewers.
- `results/`: reference outputs and supporting documents.

## Demo

The main demonstration notebook is [demo/notebook.ipynb](/Users/xiaochancai/Documents/project/github/demo/notebook.ipynb). It shows the end-to-end workflow on example systems, including:

- generation of symmetry-aware neighbor information,
- magnon-band computation,
- dispersion plotting,
- symmetry extraction,
- irreducible-representation decomposition at selected `k` points.

The `demo/` directory is intentionally tracked in Git so others can inspect example inputs and outputs without rebuilding the full workflow first.

## Example Materials

- `demo/Cu3TeO6/primitive-cell/`: example workflow for Cu3TeO6.
- `demo/RuO2/`: example workflow for RuO2.

## Current Status

This repository is still being cleaned up for public release. The core code and demo pipeline are already available locally, while documentation and notebook presentation are being refined before upload to GitHub.

## Notes

Temporary notebook checkpoints and build artifacts such as `*.o`, `*.mod`, `*.exe`, and `magnon-solver/compile/` are excluded from version control.
