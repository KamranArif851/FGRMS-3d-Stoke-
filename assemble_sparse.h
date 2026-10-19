static void assemble_iter(Model *m) {
  int nx = m->c.n[0], ny = m->c.n[1], nz = m->c.n[2], order = m->c.order,
      cell[3], nvloc[MAXV], nploc[MAXP];
  int size = m->nf + m->np + 1;
  double q[QMAX], w[QMAX];
  clock_t start = clock();
  if (size > m->c.max_equations)
    fail("Equation safety limit exceeded; reduce mesh or increase "
         "max_equations within memory capacity");
  m->neq = size;
  make_pattern(m);
  m->rhs = allocate(size, sizeof(double));
  m->mean_weights = allocate(m->np, sizeof(double));
  gauss(order, q, w);
  for (cell[0] = 0; cell[0] < nx; cell[0]++)
    for (cell[1] = 0; cell[1] < ny; cell[1]++)
      for (cell[2] = 0; cell[2] < nz; cell[2]++) {
        double *A =
            allocate((size_t)m->local_dof * m->local_dof, sizeof(double));
        double *B = allocate((size_t)m->local_p * m->local_dof, sizeof(double));
        double F[MAXDOF] = {0}, C[MAXP] = {0};
        double MP[MAXP * MAXP] = {0};
        int a, b, i, j, p, ix, iy, iz;
        connectivity(m, cell, nvloc, nploc);
        for (ix = 0; ix < order; ix++)
          for (iy = 0; iy < order; iy++)
            for (iz = 0; iz < order; iz++) {
              double xi[3] = {q[ix], q[iy], q[iz]}, x[3], N[MAXV], D[MAXV][3],
                     M[MAXP], f[3], eta;
              double W =
                  w[ix] * w[iy] * w[iz] * m->h[0] * m->h[1] * m->h[2] / 8;
              for (i = 0; i < 3; i++)
                x[i] = (cell[i] + (xi[i] + 1) / 2) * m->h[i];
              basis(m, xi, m->h, N, D, M);
              eta = coefficients(m, cell, x, M, f);
              for (a = 0; a < m->local_v; a++) {
                for (i = 0; i < 3; i++)
                  F[3 * a + i] += W * N[a] * f[i];
                for (b = 0; b < m->local_v; b++) {
                  double dot = 0;
                  for (i = 0; i < 3; i++)
                    dot += D[a][i] * D[b][i];
                  for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                      AT(A, 3 * a + i, 3 * b + j, m->local_dof) +=
                          W * eta *
                          ((i == j ? dot : 0) + D[a][j] * D[b][i] +
                           m->c.gamma * D[a][i] * D[b][j]);
                }
                for (p = 0; p < m->local_p; p++)
                  for (i = 0; i < 3; i++)
                    AT(B, p, 3 * a + i, m->local_dof) -= W * M[p] * D[a][i];
              }
              for (p = 0; p < m->local_p; p++)
                C[p] += W * M[p];
              for (int pa = 0; pa < m->local_p; pa++)
                for (int pb = 0; pb < m->local_p; pb++)
                  MP[pa * m->local_p + pb] += W * M[pa] * M[pb] / eta;
            }
        traction(m, cell, q, w, F);
        /* Scatter and add local entries. Known boundary columns go to rhs. */
        for (a = 0; a < m->local_dof; a++) {
          int ga = 3 * nvloc[a / 3] + a % 3, row = m->map[ga];
          if (row >= 0) {
            m->rhs[row] += F[a];
            for (b = 0; b < m->local_dof; b++) {
              int gb = 3 * nvloc[b / 3] + b % 3, col = m->map[gb];
              if (col >= 0)
                (*entry(&system_matrix, row, col)) += AT(A, a, b, m->local_dof);
              else
                m->rhs[row] -= AT(A, a, b, m->local_dof) * m->ub[gb];
            }
          }
          for (p = 0; p < m->local_p; p++) {
            int pr = m->nf + nploc[p];
            if (row >= 0) {
              (*entry(&system_matrix, pr, row)) += AT(B, p, a, m->local_dof);
              (*entry(&system_matrix, row, pr)) += AT(B, p, a, m->local_dof);
            } else
              m->rhs[pr] -= AT(B, p, a, m->local_dof) * m->ub[ga];
          }
        }
        for (p = 0; p < m->local_p; p++)
          m->mean_weights[nploc[p]] += C[p];
        for (int pa = 0; pa < m->local_p; pa++)
          for (int pb = 0; pb < m->local_p; pb++)
            (*entry(&pressure_mass, nploc[pa], nploc[pb])) +=
                MP[pa * m->local_p + pb];
        free(A);
        free(B);
      }
  /* A pressure gauge is needed exactly when B^T * ones is zero. */
  m->gauge = 1;
  for (int i = 0; i < m->nf; i++) {
    double sum = 0;
    for (size_t k = system_matrix.row[i]; k < system_matrix.row[i + 1]; k++)
      if (system_matrix.col[k] >= m->nf && system_matrix.col[k] < m->nf + m->np)
        sum += system_matrix.v[k];
    if (fabs(sum) > 1e-10)
      m->gauge = 0;
  }
  if (m->gauge) {
    double flux = 0, volume = 0;
    for (int p = 0; p < m->np; p++) {
      flux += m->rhs[m->nf + p];
      volume += m->mean_weights[p];
      (*entry(&system_matrix, m->nf + p, size - 1)) = m->mean_weights[p];
      (*entry(&system_matrix, size - 1, m->nf + p)) = m->mean_weights[p];
    }
    if (fabs(flux) > 1e-9 * fmax(1, norm(m->rhs + m->nf, m->np)))
      fail("Closed boundary has incompatible total inflow/outflow");
    m->rhs[size - 1] = m->c.pressure_mean / m->P0 * volume;
  } else {
    fail("This iterative implementation requires a closed domain with a "
         "pressure gauge");
  }
  for (int i = 0; i < size; i++)
    for (size_t k = system_matrix.row[i]; k < system_matrix.row[i + 1]; k++)
      m->symmetry =
          fmax(m->symmetry, fabs(system_matrix.v[k] -
                                 get(&system_matrix, system_matrix.col[k], i)));
  m->sparse_storage_MiB =
      ((double)(system_matrix.row[m->neq] + pressure_mass.row[m->np]) *
           (sizeof(int) + sizeof(double)) +
       (double)(m->neq + m->np + 2) * sizeof(size_t)) /
      1048576;
  m->assembly_time = (double)(clock() - start) / CLOCKS_PER_SEC;
}
