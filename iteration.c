/* Generalized sparse Stokes FEM: FGMRES/Schur/AMG and analytical verification.
 */
#define main reference_solver_entry
#include "stokes_fem.c"
#undef main
#include "iterative.h"
#include "assemble_sparse.h"
static void release_pc(BlockPC *p) {
  for (int c = 0; c < 3; c++) {
    amg_free(p->velocity[c]);
    free(p->index[c]);
  }
  free(p->component);
  free(p->local);
  free(p->mc);
  free(p->mp);
  free(p->temporary);
  free(p->rhs);
  free(p->answer);
  destroy(&pressure_mass);
}
int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: iteration.exe input.cfg output_prefix\n");
    return 1;
  }
  Model m = {0};
  m.c = defaults();
  read_config(&m.c, argv[1]);
  if (!has_exact(&m))
    fail("Analytical comparison requires manufactured mode and exact or noslip "
         "walls");
  double tolerance = m.c.tolerance;
  int restart = m.c.restart, maxit = m.c.max_iterations;
  printf("Domain: %.6g x %.6g x %.6g km; elements %d x %d x %d; degree %d\n",
         m.c.length[0], m.c.length[1], m.c.length[2], m.c.n[0], m.c.n[1],
         m.c.n[2], m.c.degree);
  m.L0 = fmax(m.c.length[0], fmax(m.c.length[1], m.c.length[2])) * 1000;
  m.U0 = m.c.velocity_scale * .01 / YEAR;
  m.P0 = m.c.eta0 * m.U0 / m.L0;
  m.F0 = m.P0 / m.L0;
  m.nv = m.np = 1;
  for (int d = 0; d < 3; d++) {
    m.ell[d] = m.c.length[d] * 1000 / m.L0;
    m.h[d] = m.ell[d] / m.c.n[d];
    m.vn[d] = m.c.degree * m.c.n[d] + 1;
    m.pn[d] = (m.c.degree - 1) * m.c.n[d] + 1;
    m.nv *= m.vn[d];
    m.np *= m.pn[d];
  }
  size_t estimated =
      (size_t)3 * (m.vn[0] - 2) * (m.vn[1] - 2) * (m.vn[2] - 2) + m.np + 1;
  if (estimated > (size_t)m.c.max_equations)
    fail("Mesh exceeds max_equations; reduce nodes/elements or raise the "
         "safety limit within available memory");
  printf("Velocity nodes: %d (%d x %d x %d); pressure nodes: %d. Estimated "
         "equations: %lu\n",
         m.nv, m.vn[0], m.vn[1], m.vn[2], m.np, (unsigned long)estimated);
  m.local_v = (m.c.degree + 1) * (m.c.degree + 1) * (m.c.degree + 1);
  m.local_p = m.c.degree * m.c.degree * m.c.degree;
  m.local_dof = 3 * m.local_v;
  m.map = allocate(3 * m.nv, sizeof(int));
  m.ub = allocate(3 * m.nv, sizeof(double));
  boundaries(&m);
  char itprefix[2048];
  snprintf(itprefix, sizeof(itprefix), "%s_iter", argv[2]);

  FILE *status=output_file(argv[2],"_solver.json");
  fputs("{\"converged\":false}\n",status);fclose(status);
  assemble_iter(&m);
  printf("Sparse matrix: %d equations, %lu stored entries\n", m.neq,
         (unsigned long)system_matrix.row[m.neq]);
  fflush(stdout);
  clock_t start = clock();
  BlockPC pc = build_pc(&m, itprefix);
  double setup = (double)(clock() - start) / CLOCKS_PER_SEC;
  start = clock();
  int iterations = fgmres(&m, &pc, itprefix, tolerance, restart, maxit);
  m.solve_time = (double)(clock() - start) / CLOCKS_PER_SEC;
  /* Save the actual input so later edits cannot change interpretation of
   * results. */
  FILE *source = fopen(argv[1], "rb"),
       *snapshot = output_file(argv[2], "_input.cfg");
  if (!source)
    fail("Cannot snapshot input");
  int ch;
  while ((ch = fgetc(source)) != EOF)
    fputc(ch, snapshot);
  fclose(source);
  fclose(snapshot);
  errors_and_stats(&m, itprefix);
  export_results(&m, itprefix);
  export_all_fields(&m, itprefix);
  double iterative_time = m.solve_time;
  release_pc(&pc);
  FILE *f = output_file(argv[2], "_solver.json");
  double *res = allocate(m.neq, sizeof(double));
  multiply(&system_matrix, m.solution, res);
  for (int i = 0; i < m.neq; i++)
    res[i] = m.rhs[i] - res[i];
  double true_res = norm(res, m.neq) / fmax(norm(m.rhs, m.neq), 1e-30);
  free(res);
  fprintf(
      f,
      "{\n\"converged\":true,\n\"method\":\"FGMRES\",\n\"velocity_nodes\":%d,\n\"pressure_nodes\":%"
      "d,\n\"equations\":%d,\n\"iterations\":%d,\n\"relative_tolerance\":%.17g,"
      "\n\"restart\":%d,\n\"max_iterations\":%d,\n\"true_relative_residual\":%."
      "17g,\n\"setup_seconds\":%.17g,\n\"solve_seconds\":%.17g\n}\n",
      m.nv, m.np, m.neq, iterations, tolerance, restart, maxit, true_res, setup,
      iterative_time);
  fclose(f);
  printf("FGMRES converged in %d iterations. True relative residual %.3e. "
         "Setup %.3fs; solve %.3fs.\n",
         iterations, true_res, setup, iterative_time);
  destroy(&system_matrix);
  free(m.solution);
  free(m.velocity);
  free(m.pressure);
  free(m.map);
  free(m.ub);
  free(m.rhs);
  free(m.mean_weights);
  return 0;
}
