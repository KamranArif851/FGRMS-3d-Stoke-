/* Genuine algebraic aggregation multigrid, CSR storage and flexible GMRES. */
typedef struct {
  int n;
  size_t *row;
  int *col;
  double *v;
} Sparse;
typedef struct {
  int *a, n, cap;
} PatternRow;
static Sparse system_matrix, pressure_mass;

static void push(PatternRow *r, int col) {
  if (r->n == r->cap) {
    r->cap = r->cap ? 2 * r->cap : 32;
    int *p = realloc(r->a, (size_t)r->cap * sizeof(int));
    if (!p)
      fail("Pattern allocation");
    r->a = p;
  }
  r->a[r->n++] = col;
}
static int compare_int(const void *a, const void *b) {
  int x = *(const int *)a, y = *(const int *)b;
  return (x > y) - (x < y);
}
static Sparse compress_pattern(PatternRow *rows, int n) {
  Sparse a = {0};
  a.n = n;
  a.row = allocate(n + 1, sizeof(size_t));
  for (int i = 0; i < n; i++) {
    qsort(rows[i].a, rows[i].n, sizeof(int), compare_int);
    int unique = 0;
    for (int k = 0; k < rows[i].n; k++)
      if (!k || rows[i].a[k] != rows[i].a[k - 1])
        rows[i].a[unique++] = rows[i].a[k];
    rows[i].n = unique;
    a.row[i + 1] = a.row[i] + unique;
  }
  a.col = allocate(a.row[n], sizeof(int));
  a.v = allocate(a.row[n], sizeof(double));
  for (int i = 0; i < n; i++) {
    memcpy(a.col + a.row[i], rows[i].a, rows[i].n * sizeof(int));
    free(rows[i].a);
  }
  free(rows);
  return a;
}
static double *entry(Sparse *a, int i, int j) {
  size_t lo = a->row[i], hi = a->row[i + 1];
  while (lo < hi) {
    size_t m = (lo + hi) / 2;
    if (a->col[m] < j)
      lo = m + 1;
    else
      hi = m;
  }
  if (lo == a->row[i + 1] || a->col[lo] != j)
    fail("Sparse pattern missing entry");
  return a->v + lo;
}
static double get(const Sparse *a, int i, int j) {
  size_t lo = a->row[i], hi = a->row[i + 1];
  while (lo < hi) {
    size_t m = (lo + hi) / 2;
    if (a->col[m] < j)
      lo = m + 1;
    else
      hi = m;
  }
  return lo < a->row[i + 1] && a->col[lo] == j ? a->v[lo] : 0;
}
static void multiply(const Sparse *a, const double *x, double *y) {
  for (int i = 0; i < a->n; i++) {
    double sum = 0;
    for (size_t k = a->row[i]; k < a->row[i + 1]; k++)
      sum += a->v[k] * x[a->col[k]];
    y[i] = sum;
  }
}
static void destroy(Sparse *a) {
  free(a->row);
  free(a->col);
  free(a->v);
  memset(a, 0, sizeof(*a));
}
static double dot_product(const double *a, const double *b, int n) {
  double s = 0;
  for (int i = 0; i < n; i++)
    s += a[i] * b[i];
  return s;
}
static void make_pattern(Model *m) {
  int size = m->nf + m->np + 1;
  PatternRow *r = allocate(size, sizeof(*r));
  PatternRow *mass = allocate(m->np, sizeof(*mass));
  int cell[3], v[MAXV], p[MAXP];
  for (cell[0] = 0; cell[0] < m->c.n[0]; cell[0]++)
    for (cell[1] = 0; cell[1] < m->c.n[1]; cell[1]++)
      for (cell[2] = 0; cell[2] < m->c.n[2]; cell[2]++) {
        connectivity(m, cell, v, p);
        for (int a = 0; a < m->local_p; a++)
          for (int b = 0; b < m->local_p; b++)
            push(mass + p[a], p[b]);
        int list[MAXDOF], count = 0;
        for (int a = 0; a < m->local_dof; a++) {
          int row = m->map[3 * v[a / 3] + a % 3];
          if (row >= 0)
            list[count++] = row;
        }
        for (int a = 0; a < count; a++) {
          for (int b = 0; b < count; b++)
            push(r + list[a], list[b]);
          for (int b = 0; b < m->local_p; b++) {
            push(r + list[a], m->nf + p[b]);
            push(r + m->nf + p[b], list[a]);
          }
        }
      }
  for (int p = 0; p < m->np; p++) {
    push(r + m->nf + p, size - 1);
    push(r + size - 1, m->nf + p);
  }
  system_matrix = compress_pattern(r, size);
  pressure_mass = compress_pattern(mass, m->np);
}
typedef struct AMG {
  Sparse a;
  int *aggregate, nc, level;
  double *diagonal, *scratch, *coarse_rhs, *coarse_x, *chol;
  struct AMG *coarse;
} AMG;
static void cholesky(double *a, int n) {
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
      double s = AT(a, i, j, n);
      for (int k = 0; k < j; k++)
        s -= AT(a, i, k, n) * AT(a, j, k, n);
      if (i == j) {
        if (s <= 0)
          fail("Nonpositive coarse matrix");
        AT(a, i, j, n) = sqrt(s);
      } else
        AT(a, i, j, n) = s / AT(a, j, j, n);
    }
}
static void chol_solve(const double *a, int n, const double *b, double *x) {
  for (int i = 0; i < n; i++) {
    double s = b[i];
    for (int j = 0; j < i; j++)
      s -= AT(a, i, j, n) * x[j];
    x[i] = s / AT(a, i, i, n);
  }
  for (int i = n - 1; i >= 0; i--) {
    double s = x[i];
    for (int j = i + 1; j < n; j++)
      s -= AT(a, j, i, n) * x[j];
    x[i] = s / AT(a, i, i, n);
  }
}
static AMG *amg_build(Sparse a, int level, FILE *log) {
  AMG *g = allocate(1, sizeof(*g));
  g->a = a;
  g->level = level;
  int n = a.n;
  g->diagonal = allocate(n, sizeof(double));
  g->scratch = allocate(n, sizeof(double));
  for (int i = 0; i < n; i++) {
    g->diagonal[i] = get(&a, i, i);
    if (g->diagonal[i] <= 0)
      fail("AMG requires positive diagonal");
  }
  fprintf(log, "%d,%d,%lu\n", level, n, (unsigned long)a.row[n]);
  if (n <= 32 || level >= 12) {
    g->chol = allocate((size_t)n * n, sizeof(double));
    for (int i = 0; i < n; i++)
      for (size_t k = a.row[i]; k < a.row[i + 1]; k++)
        AT(g->chol, i, a.col[k], n) = a.v[k];
    cholesky(g->chol, n);
    return g;
  }
  g->aggregate = allocate(n, sizeof(int));
  for (int i = 0; i < n; i++)
    g->aggregate[i] = -1;
  int nc = 0;
  for (int i = 0; i < n; i++)
    if (g->aggregate[i] < 0) {
      int partner = -1;
      double strongest = 0;
      for (size_t k = a.row[i]; k < a.row[i + 1]; k++) {
        int j = a.col[k];
        double score = fabs(a.v[k]) / sqrt(g->diagonal[i] * g->diagonal[j]);
        if (j != i && g->aggregate[j] < 0 && score > strongest) {
          strongest = score;
          partner = j;
        }
      }
      g->aggregate[i] = nc;
      if (partner >= 0)
        g->aggregate[partner] = nc;
      nc++;
    }
  g->nc = nc;
  PatternRow *rows = allocate(nc, sizeof(*rows));
  for (int i = 0; i < n; i++)
    for (size_t k = a.row[i]; k < a.row[i + 1]; k++)
      push(rows + g->aggregate[i], g->aggregate[a.col[k]]);
  Sparse ac = compress_pattern(rows, nc);
  for (int i = 0; i < n; i++)
    for (size_t k = a.row[i]; k < a.row[i + 1]; k++)
      *entry(&ac, g->aggregate[i], g->aggregate[a.col[k]]) += a.v[k];
  g->coarse_rhs = allocate(nc, sizeof(double));
  g->coarse_x = allocate(nc, sizeof(double));
  g->coarse = amg_build(ac, level + 1, log);
  return g;
}
static void sweep(AMG *g, const double *b, double *x, int reverse) {
  int n = g->a.n;
  for (int t = 0; t < n; t++) {
    int i = reverse ? n - 1 - t : t;
    double s = b[i];
    for (size_t k = g->a.row[i]; k < g->a.row[i + 1]; k++)
      if (g->a.col[k] != i)
        s -= g->a.v[k] * x[g->a.col[k]];
    x[i] = s / g->diagonal[i];
  }
}
static void cycle(AMG *g, const double *b, double *x) {
  if (g->chol) {
    chol_solve(g->chol, g->a.n, b, x);
    return;
  }
  for (int s = 0; s < 2; s++)
    sweep(g, b, x, 0);
  multiply(&g->a, x, g->scratch);
  memset(g->coarse_rhs, 0, g->nc * sizeof(double));
  memset(g->coarse_x, 0, g->nc * sizeof(double));
  for (int i = 0; i < g->a.n; i++)
    g->coarse_rhs[g->aggregate[i]] += b[i] - g->scratch[i];
  cycle(g->coarse, g->coarse_rhs, g->coarse_x);
  for (int i = 0; i < g->a.n; i++)
    x[i] += g->coarse_x[g->aggregate[i]];
  for (int s = 0; s < 2; s++)
    sweep(g, b, x, 1);
}
static void amg_free(AMG *g) {
  if (!g)
    return;
  amg_free(g->coarse);
  destroy(&g->a);
  free(g->aggregate);
  free(g->diagonal);
  free(g->scratch);
  free(g->coarse_rhs);
  free(g->coarse_x);
  free(g->chol);
  free(g);
}
/* Inner preconditioned CG is applied only to SPD subblocks, never to Stokes K.
 */
static void pcg(const Sparse *a, AMG *g, const double *b, double *x, double tol,
                int maxit) {
  int n = a->n;
  double *r = allocate(n, sizeof(double)), *z = allocate(n, sizeof(double)),
         *p = allocate(n, sizeof(double)), *ap = allocate(n, sizeof(double));
  memset(x, 0, n * sizeof(double));
  memcpy(r, b, n * sizeof(double));
  double target = tol * norm(b, n), old = 0;
  for (int it = 0; it < maxit && norm(r, n) > fmax(target, 1e-30); it++) {
    memset(z, 0, n * sizeof(double));
    if (g)
      cycle(g, r, z);
    else
      for (int i = 0; i < n; i++)
        z[i] = r[i] / get(a, i, i);
    double rz = dot_product(r, z, n);
    if (it == 0)
      memcpy(p, z, n * sizeof(double));
    else
      for (int i = 0; i < n; i++)
        p[i] = z[i] + rz / old * p[i];
    multiply(a, p, ap);
    double pap = dot_product(p, ap, n);
    if (pap <= 0)
      fail("Nonpositive inner CG curvature");
    double alpha = rz / pap;
    for (int i = 0; i < n; i++) {
      x[i] += alpha * p[i];
      r[i] -= alpha * ap[i];
    }
    old = rz;
  }
  free(r);
  free(z);
  free(p);
  free(ap);
}
typedef struct {
  Model *m;
  AMG *velocity[3];
  int *index[3], *component, *local, count[3];
  double *mc, *mp, *temporary, *rhs, *answer;
  double cMc;
} BlockPC;
static BlockPC build_pc(Model *m, const char *prefix) {
  BlockPC p = {0};
  p.m = m;
  p.component = allocate(m->nf, sizeof(int));
  p.local = allocate(m->nf, sizeof(int));
  for (int i = 0; i < 3 * m->nv; i++)
    if (m->map[i] >= 0) {
      int r = m->map[i], c = i % 3;
      p.component[r] = c;
      p.local[r] = p.count[c]++;
    }
  for (int c = 0; c < 3; c++) {
    p.index[c] = allocate(p.count[c], sizeof(int));
    for (int i = 0; i < m->nf; i++)
      if (p.component[i] == c)
        p.index[c][p.local[i]] = i;
    PatternRow *r = allocate(p.count[c], sizeof(*r));
    for (int i = 0; i < p.count[c]; i++) {
      int gi = p.index[c][i];
      for (size_t k = system_matrix.row[gi]; k < system_matrix.row[gi + 1];
           k++) {
        int j = system_matrix.col[k];
        if (j < m->nf && p.component[j] == c)
          push(r + i, p.local[j]);
      }
    }
    Sparse a = compress_pattern(r, p.count[c]);
    for (int i = 0; i < a.n; i++)
      for (size_t k = a.row[i]; k < a.row[i + 1]; k++)
        a.v[k] = get(&system_matrix, p.index[c][i], p.index[c][a.col[k]]);
    char suffix[64];
    snprintf(suffix, sizeof(suffix), "_amg_%c.csv", 'x' + c);
    FILE *log = output_file(prefix, suffix);
    fprintf(log, "level,unknowns,nonzeros\n");
    p.velocity[c] = amg_build(a, 0, log);
    fclose(log);
  }

  p.mc = allocate(m->np, sizeof(double));
  p.mp = allocate(m->np, sizeof(double));
  p.temporary = allocate(m->nf, sizeof(double));
  p.rhs = allocate(m->nf, sizeof(double));
  p.answer = allocate(m->nf, sizeof(double));
  pcg(&pressure_mass, NULL, m->mean_weights, p.mc, 1e-12, 200);
  p.cMc = dot_product(m->mean_weights, p.mc, m->np);
  if (!isfinite(p.cMc) || p.cMc <= 0)
    fail("Invalid pressure gauge preconditioner");
  return p;
}
static void precondition(BlockPC *p, const double *r, double *z) {
  Model *m = p->m;
  pcg(&pressure_mass, NULL, r + m->nf, p->mp, m->c.pressure_tolerance,
      m->c.pressure_iterations);
  double lambda =
      (r[m->neq - 1] + dot_product(m->mean_weights, p->mp, m->np)) / p->cMc;
  z[m->neq - 1] = lambda;
  for (int i = 0; i < m->np; i++)
    z[m->nf + i] = -p->mp[i] + p->mc[i] * lambda;
  for (int i = 0; i < m->nf; i++) {
    double b = r[i];
    for (size_t k = system_matrix.row[i]; k < system_matrix.row[i + 1]; k++) {
      int j = system_matrix.col[k];
      if (j >= m->nf && j < m->nf + m->np)
        b -= system_matrix.v[k] * z[j];
    }
    p->temporary[i] = b;
  }
  for (int c = 0; c < 3; c++) {
    for (int i = 0; i < p->count[c]; i++)
      p->rhs[i] = p->temporary[p->index[c][i]];
    pcg(&p->velocity[c]->a, p->velocity[c], p->rhs, p->answer,
        m->c.velocity_tolerance, m->c.velocity_iterations);
    for (int i = 0; i < p->count[c]; i++)
      z[p->index[c][i]] = p->answer[i];
  }
}
static int fgmres(Model *m, BlockPC *pc, const char *prefix, double tolerance,
                  int restart, int maxit) {
  int n = m->neq, total = 0;
  double *x = allocate(n, sizeof(double)),
         *v = allocate((size_t)(restart + 1) * n, sizeof(double)),
         *z = allocate((size_t)restart * n, sizeof(double)),
         *h = allocate((size_t)(restart + 1) * restart, sizeof(double)),
         *cs = allocate(restart, sizeof(double)),
         *sn = allocate(restart, sizeof(double)),
         *g = allocate(restart + 1, sizeof(double)),
         *y = allocate(restart, sizeof(double)),
         *r = allocate(n, sizeof(double)), *w = allocate(n, sizeof(double));
  double bn = norm(m->rhs, n), den = fmax(bn, 1e-30);
  FILE *log = output_file(prefix, "_history.csv");
  fprintf(log,
          "iteration,estimated_relative_residual,true_relative_residual\n");
  while (total < maxit) {
    multiply(&system_matrix, x, w);
    for (int i = 0; i < n; i++)
      r[i] = m->rhs[i] - w[i];
    double beta = norm(r, n);
    fprintf(log, "%d,%.17g,%.17g\n", total, beta / den, beta / den);
    if (beta <= tolerance * den)
      break;
    memset(h, 0, (size_t)(restart + 1) * restart * sizeof(double));
    memset(g, 0, (restart + 1) * sizeof(double));
    g[0] = beta;
    for (int i = 0; i < n; i++)
      v[i] = r[i] / beta;
    int used = 0;
    for (int j = 0; j < restart && total < maxit; j++) {
      precondition(pc, v + (size_t)j * n, z + (size_t)j * n);
      multiply(&system_matrix, z + (size_t)j * n, w);
      /* Twice modified Gram-Schmidt limits loss of orthogonality. */
      for (int pass = 0; pass < 2; pass++)
        for (int i = 0; i <= j; i++) {
          double a = dot_product(w, v + (size_t)i * n, n);
          AT(h, i, j, restart) += a;
          for (int k = 0; k < n; k++)
            w[k] -= a * v[(size_t)i * n + k];
        }
      double next = norm(w, n);
      AT(h, j + 1, j, restart) = next;
      if (next > 1e-30)
        for (int k = 0; k < n; k++)
          v[(size_t)(j + 1) * n + k] = w[k] / next;
      for (int i = 0; i < j; i++) {
        double a = AT(h, i, j, restart), b = AT(h, i + 1, j, restart);
        AT(h, i, j, restart) = cs[i] * a + sn[i] * b;
        AT(h, i + 1, j, restart) = -sn[i] * a + cs[i] * b;
      }
      double a = AT(h, j, j, restart), b = AT(h, j + 1, j, restart),
             rho = hypot(a, b);
      if (rho < 1e-30)
        fail("FGMRES breakdown");
      cs[j] = a / rho;
      sn[j] = b / rho;
      AT(h, j, j, restart) = rho;
      AT(h, j + 1, j, restart) = 0;
      g[j + 1] = -sn[j] * g[j];
      g[j] *= cs[j];
      total++;
      used = j + 1;
      fprintf(log, "%d,%.17g,nan\n", total, fabs(g[j + 1]) / den);
      if (total % 20 == 0) {
        printf("FGMRES %d estimated relative residual %.3e\n", total,
               fabs(g[j + 1]) / den);
        fflush(stdout);
      }
      if (fabs(g[j + 1]) <= tolerance * den || next < 1e-30)
        break;
    }
    for (int i = used - 1; i >= 0; i--) {
      double s = g[i];
      for (int j = i + 1; j < used; j++)
        s -= AT(h, i, j, restart) * y[j];
      y[i] = s / AT(h, i, i, restart);
    }
    for (int j = 0; j < used; j++)
      for (int i = 0; i < n; i++)
        x[i] += z[(size_t)j * n + i] * y[j];
  }
  multiply(&system_matrix, x, w);
  for (int i = 0; i < n; i++)
    r[i] = m->rhs[i] - w[i];
  double true_rel = norm(r, n) / den;
  fprintf(log, "%d,%.17g,%.17g\n", total, true_rel, true_rel);
  fclose(log);
  if (!isfinite(true_rel) || true_rel > tolerance)
    fail("FGMRES failed true residual tolerance");
  m->solution = x;
  m->residual = norm(r, n) / fmax(1, bn);
  m->lambda = x[n - 1];
  double cont = 0;
  for (int i = 0; i < m->np; i++) {
    double c = r[m->nf + i] + m->mean_weights[i] * m->lambda;
    cont += c * c;
  }
  m->continuity = sqrt(cont);
  m->velocity = allocate(3 * m->nv, sizeof(double));
  m->pressure = allocate(m->np, sizeof(double));
  for (int i = 0; i < 3 * m->nv; i++)
    m->velocity[i] = m->map[i] >= 0 ? x[m->map[i]] : m->ub[i];
  memcpy(m->pressure, x + m->nf, m->np * sizeof(double));
  free(v);
  free(z);
  free(h);
  free(cs);
  free(sn);
  free(g);
  free(y);
  free(r);
  free(w);
  return total;
}
