#Generalized 3-D Stokes Verification in C

Overview

This program implements a generalized three-dimensional Stokes verification framework in C to evaluate the accuracy and convergence of finite element solutions against analytical benchmarks.

To run the program, review START_HERE.txt, configure the parameters in input.cfg, and execute run.bat. The generated results can be visualized using gallery.html, which compares analytical solutions (blue) with numerical solutions (red).

Numerical Methodology

The Stokes system is solved using restarted Flexible Generalized Minimal Residual (FGMRES) with a block Schur-complement preconditioner. Algebraic aggregation-based multigrid (AMG) accelerates the velocity subproblems, while an inverse-viscosity-weighted pressure mass matrix approximates the pressure Schur complement. The global Stokes and pressure mass matrices are assembled directly in compressed sparse row (CSR) format.

The framework supports polynomial degrees (p=2,3,4) and configurable domain dimensions, mesh resolution, physical parameters, integration orders, and solver tolerances through input.cfg.

Scope and Outputs

The program is designed for manufactured-solution verification on rectangular three-dimensional domains under supported boundary conditions and viscosity distributions. It is intended for numerical validation rather than unrestricted geodynamic simulation.

Key outputs include:

results/main_solver.json: Solver convergence and performance diagnostics.

results/main_iter_stats.json: Finite element errors and numerical checks.

results/main_iter_history.csv: Iterative convergence history.

results/main_iter_amg_*.csv: AMG hierarchy information.

Additional CSV and VTK files: Numerical fields, stresses, and profiles.

figures/ and gallery.html: Visual comparisons of analytical and numerical solutions.

study_gallery.html and node_comparison.csv: Mesh-refinement comparison results.

Reported error measures are dimensionless unless physical units are explicitly specified. The true relative residual is defined as

[
r_{\mathrm{true}}=
\frac{|\mathbf b-K\mathbf x|_2}
{\max(|\mathbf b|_2,10^{-30})}.
]

The reported sparse matrix storage excludes AMG hierarchies, temporary structures, and other solver allocations, and therefore does not represent total peak memory usage.

Code Organization

iteration.c: Input processing, assembly, iterative solution, and output.

stokes_fem.c: Finite element basis functions, analytical solutions, boundary conditions, and error evaluation.

assemble_sparse.h: Sparse finite element assembly.

iterative.h: Sparse algebra, AMG, inner CG, block preconditioning, and FGMRES.

plots.c: Field evaluation and SVG plotting.

comparison_plots.c: Solution comparisons, residual plots, and mesh studies.

verify_coefficients.c: Optional regression checking.

Benchmark Validation

The reference benchmark uses a (1000\ \mathrm{km}) cubic domain, a Q4/Q3 finite element pair, and a (4\times4\times3) element mesh. The discretization contains 3,757 velocity nodes, 1,690 pressure nodes, and 9,116 algebraic equations.

The documented results show convergence in 69 FGMRES iterations, with a true relative residual of approximately (9.32\times10^{-13}). The relative (L_2) errors are approximately 0.0442455% for velocity and 0.150469% for pressure.

These values provide reference results for validating the implementation under the same mesh, physical parameters, boundary conditions, and solver settings. Additional verification can be performed through mesh refinement, residual analysis, and comparison of numerical and analytical fields.

For further details, consult VALIDATION.txt. Solver performance may vary with mesh resolution, domain aspect ratio, viscosity distribution, and convergence settings.
