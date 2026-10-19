# Generalized 3-D Stokes verification in C

Read **START_HERE.txt**, edit **input.cfg**, then double-click **run.bat**.
Open **gallery.html** for analytical (blue) versus iterative (red) plots.

This version removes the global LU solver, dense global assembly and LU comparisons. It uses restarted right FGMRES with a block Schur preconditioner, algebraic aggregation multigrid for the component velocity blocks, and an inverse-viscosity weighted pressure mass approximation. Global Stokes and pressure mass matrices are both assembled directly in CSR. Small coarse AMG systems use Cholesky; this is part of AMG, not a comparison solve of the full Stokes system.

FGMRES is the flexible variant of GMRES. Approximate inner CG solves depend on the current right-hand side and can use different numbers of iterations, so ordinary fixed-preconditioner GMRES is not substituted silently.

## Editable inputs

All numerical settings are in input.cfg. Domain lengths, mesh description, polynomial degree, physical scales, pressure mean, benchmark, integration orders, outer tolerance/restart/iteration limit and inner-block tolerances/limits can be changed there. The chosen mesh description is authoritative: `mesh_input = elements` uses `elements`; `mesh_input = nodes` uses `velocity_nodes` and derives the element counts.

For degree p:

- Velocity nodes: (p*nx+1)(p*ny+1)(p*nz+1).
- Pressure nodes: ((p-1)*nx+1)((p-1)*ny+1)((p-1)*nz+1).
- Each axis must have at least two elements; p is 2, 3 or 4.

This is a manufactured-solution verification program on rectangular boxes, not an unrestricted geodynamic simulator. Exact/noslip walls and constant/linear unclipped viscosity preserve the reference test. Different physical boundary data require a new matching analytical solution and forcing. Such unsupported comparisons are rejected before solving.

## Outputs

- results/main_input.cfg: exact input snapshot used by the solver.
- results/main_solver.json: actual outer iterations, requested settings, true relative residual, setup/solve timings.
- results/main_iter_stats.json: integrated FEM errors, divergence and other checks.
- results/main_iter_history.csv: Arnoldi estimates and true residual checks. `nan` means the true residual was not recomputed at that intermediate step.
- results/main_iter_amg_x.csv (and y,z): hierarchy levels.
- Other main_iter CSV and VTK files: nodal fields, sampled fields, stresses and profiles.
- figures/: SVG images; gallery.html: single-case results.
- study_gallery.html: delivered three-mesh gallery; node_comparison.csv: mesh error table.

Errors in stats.json are dimensionless unless a key explicitly states physical units. Plot labels give physical units. Stress norms count symmetric off-diagonal components twice. Three-dimensional surface heights represent field values on a z/Lz=0.43 slice, not physical vertical position. Actual 3-D flow is shown by streamlines.

The true relative stopping residual is ||b-Kx||/max(||b||,1e-30). The legacy `scaled_residual` statistic uses max(||b||,1); do not confuse these denominators. `sparse_matrix_storage_MiB` includes CSR storage for the global Stokes matrix and pressure mass matrix; it excludes AMG levels, temporary assembly patterns, Krylov vectors and other allocations, so it is not peak process memory.

## Code map

- iteration.c: read inputs, resolve mesh, assemble, precondition, solve and export.
- stokes_fem.c: shared FEM basis, exact solution, force, boundary treatment, error integration and exports. 
- assemble_sparse.h: direct CSR finite-element assembly.
- iterative.h: sparse algebra, AMG, inner CG, block preconditioner, restarted FGMRES.
- plots.c: C-only field evaluation and SVG drawing functions.
- comparison_plots.c: analytical/iterative gallery, residual histories and optional mesh study.
- verify_coefficients.c: optional regression checker for iterative CSVs; not required to run the model.

The previous Iteration Method folder is preserved and contains the earlier full mathematical PDF and historical LU comparisons. That historical PDF describes the earlier version; use this folder's instructions for current settings and outputs.

## Validation

See VALIDATION.txt for tests. Default 1000 km cube, Q4/Q3, 4x4x3 elements: 3757 velocity nodes, 1690 pressure nodes, 9116 equations. FGMRES converges in 69 iterations to about 9.32e-13 true relative residual. Relative L2 FEM errors are approximately 0.0442455% velocity and 0.150469% pressure. Errors match the earlier iterative implementation; sparse mass assembly and removing LU do not change the benchmark solution.

All numerical computation and plotting is C. Compiler flags enable warnings as errors. No claim is made that every aspect ratio, viscosity distribution or arbitrarily large mesh will converge cheaply; the supplied component AMG is a compact verification implementation.
