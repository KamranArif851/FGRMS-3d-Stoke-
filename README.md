# Generalized 3D Stokes Verification in C

## Overview

This project implements a generalized three-dimensional Stokes solver in C to verify finite element solutions against analytical benchmarks.

The numerical framework uses restarted Flexible Generalized Minimal Residual (FGMRES) with a block Schur-complement preconditioner, algebraic aggregation-based multigrid (AMG) for velocity subproblems, and an inverse-viscosity-weighted pressure mass approximation. The global Stokes and pressure mass matrices are assembled directly in compressed sparse row (CSR) format.

## Getting Started

1. Read `START_HERE.txt`.
2. Configure the simulation parameters in `input.cfg`.
3. Run `run.bat` to compile and execute the solver.
4. Open `gallery.html` to compare analytical (blue) and numerical (red) solutions.

## Numerical Methodology

The steady Stokes equations for an incompressible viscous fluid are

$$
-\nabla \cdot \boldsymbol{\sigma} = \rho \mathbf{g}
$$

$$
\nabla \cdot \mathbf{u} = 0
$$

where:

- $\boldsymbol{\sigma}$ is the stress tensor.
- $\rho$ is the density.
- $\mathbf{g}$ is gravitational acceleration.
- $\mathbf{u}$ is the velocity field.

For a Newtonian viscous material, the stress tensor can be expressed as

$$
\boldsymbol{\sigma}
= -p\mathbf{I} + 2\eta\boldsymbol{\dot{\varepsilon}}
$$

where $p$ is pressure, $\mathbf{I}$ is the identity tensor, $\eta$ is viscosity, and $\boldsymbol{\dot{\varepsilon}}$ is the strain-rate tensor.

The solver supports polynomial degrees $p=2,3,4$ and configurable mesh dimensions, physical parameters, integration orders, and convergence tolerances.

For polynomial degree $p$, the numbers of velocity and pressure nodes are

$$
N_v=(pn_x+1)(pn_y+1)(pn_z+1)
$$

$$
N_p=((p-1)n_x+1)((p-1)n_y+1)((p-1)n_z+1)
$$

where $n_x$, $n_y$, and $n_z$ are the element counts along the three coordinate directions.

## Project Structure

| File | Description |
|---|---|
| `iteration.c` | Main solver workflow and output generation |
| `stokes_fem.c` | Finite element formulation, analytical solutions, and error evaluation |
| `assemble_sparse.h` | Sparse finite element assembly |
| `iterative.h` | FGMRES, AMG, inner CG, and block preconditioning |
| `plots.c` | Field evaluation and SVG plotting |
| `comparison_plots.c` | Analytical comparisons, residual plots, and mesh studies |
| `input.cfg` | Simulation configuration |
| `VALIDATION.txt` | Validation tests and reference results |

## Outputs

The program generates solver diagnostics, convergence histories, finite element error statistics, field data, and visualizations.

- `results/main_solver.json` — Solver iterations, true relative residual, and timings.
- `results/main_iter_stats.json` — Numerical errors and diagnostic measures.
- `results/main_iter_history.csv` — Iterative convergence history.
- `results/main_iter_amg_x.csv` — AMG hierarchy information for the x-velocity block, with corresponding files for y and z.
- `gallery.html` — Analytical and numerical solution comparisons.
- `study_gallery.html` — Three-mesh comparison gallery.
- `node_comparison.csv` — Mesh error comparison.
- `figures/` — SVG plots and visualizations.

## Benchmark Validation

The reference benchmark uses a $1000\ \mathrm{km}$ cubic domain with a Q4/Q3 finite element pair and a $4\times4\times3$ element mesh.

| Quantity | Reference result |
|---|---:|
| Velocity nodes | 3,757 |
| Pressure nodes | 1,690 |
| Algebraic equations | 9,116 |
| FGMRES iterations | 69 |
| True relative residual | $9.32\times10^{-13}$ |
| Relative velocity $L_2$ error | 0.0442455% |
| Relative pressure $L_2$ error | 0.150469% |

The true relative residual is calculated as

$$
r_{\mathrm{true}} =
\frac{\|\mathbf{b}-K\mathbf{x}\|_2}
{\max(\|\mathbf{b}\|_2,10^{-30})}
$$

where $K$ is the assembled Stokes matrix, $\mathbf{x}$ is the numerical solution, and $\mathbf{b}$ is the right-hand-side vector.

These reference values can be used to assess reproducibility under the same mesh, physical parameters, boundary conditions, and solver tolerances.

## Limitations

This implementation is intended for manufactured-solution verification on rectangular three-dimensional domains, rather than unrestricted geodynamic simulation.

The supported benchmarks require consistent analytical solutions, forcing terms, boundary conditions, and viscosity distributions. Solver convergence and computational cost may vary with mesh resolution, domain aspect ratio, viscosity distribution, and numerical settings.

See `VALIDATION.txt` for additional validation details.

## License

No license has been specified. Add a license file if you intend to distribute the project under particular reuse conditions.
