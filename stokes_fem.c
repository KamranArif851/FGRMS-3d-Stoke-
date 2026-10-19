/* Shared Qp/Q(p-1) FEM routines. Used by iteration.c and comparison_plots.c. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
static void *virtual_blocks[32];
static void release_memory(void *p) {
  for (int k = 0; k < 32; k++)
    if (p && virtual_blocks[k] == p) {
      VirtualFree(p, 0, MEM_RELEASE);
      virtual_blocks[k] = NULL;
      return;
    }
  free(p);
}
#define free release_memory
#endif

#define QMAX 10
#define MAXV 125
#define MAXP 64
#define MAXDOF 375
#define YEAR (365.25 * 24.0 * 3600.0)
#define AT(a, i, j, n) ((a)[(size_t)(i) * (n) + (j)])
enum { EXACT, NOSLIP, FREESLIP, VELOCITY, TRACTION };
typedef struct {
  int kind;
  double value[3];
} Boundary;
typedef struct {
  int n[3], order, error_order, max_equations, mode, viscosity, averaging,
      degree, benchmark, mesh_nodes, requested_nodes[3], restart,
      max_iterations, velocity_iterations, pressure_iterations;
  double tolerance, velocity_tolerance, pressure_tolerance;
  double length[3], eta0, velocity_scale, pressure_mean, gamma;
  double density, anomaly, gravity[3], center[3], width, contrast, eta_min,
      eta_max;
  Boundary bc[6];
} Config;
typedef struct {
  Config c;
  int vn[3], pn[3], nv, np, nf, neq, gauge, local_v, local_p, local_dof;
  double L0, U0, P0, F0, ell[3], h[3];
  int *map;
  double *ub, *velocity, *pressure, *rhs, *solution, *mean_weights;
  double residual, continuity, boundary_error, lambda, symmetry, assembly_time,
      solve_time, sparse_storage_MiB;
} Model;
static const char *face_names[6] = {"x_min", "x_max", "y_min",
                                    "y_max", "z_min", "z_max"};
static void fail(const char *s) {
  fprintf(stderr, "ERROR: %s\n", s);
  exit(EXIT_FAILURE);
}
static void *allocate(size_t n, size_t size) {
  void *p;
  if (size && n > SIZE_MAX / size)
    fail("Allocation size overflow");
  p = NULL;
#ifdef _WIN32
  if (n * size >= 64u * 1024u * 1024u) {
    for (int k = 0; k < 32; k++)
      if (!virtual_blocks[k]) {
        p = VirtualAlloc(NULL, n * size, MEM_RESERVE | MEM_COMMIT,
                         PAGE_READWRITE);
        virtual_blocks[k] = p;
        break;
      }
  } else
#endif
    p = calloc(n, size);
  if (!p)
    fail("Insufficient memory; reduce element counts");
  return p;
}
static FILE *output_file(const char *prefix, const char *suffix) {
  char path[2048];
  FILE *f;
  if (snprintf(path, sizeof(path), "%s%s", prefix, suffix) >= (int)sizeof(path))
    fail("Output path too long");
  f = fopen(path, "w");
  if (!f)
    fail("Cannot open output; create the output directory first");
  return f;
}
static int node_id(int i, int j, int k, const int count[3]) {
  return (i * count[1] + j) * count[2] + k;
}
static double norm(const double *a, int n) {
  double s = 0;
  int i;
  for (i = 0; i < n; i++)
    s += a[i] * a[i];
  return sqrt(s);
}
static Config defaults(void) {
  Config c;
  int d, f;
  memset(&c, 0, sizeof(c));
  for (d = 0; d < 3; d++) {
    c.n[d] = 4;
    c.length[d] = 400;
    c.center[d] = 0.5;
  }
  c.center[1] = 0.45;
  c.order = 5;
  c.error_order = 6;
  c.max_equations = 14000;
  c.degree = 2;
  c.tolerance = 1e-12;
  c.restart = 60;
  c.max_iterations = 600;
  c.velocity_tolerance = 1e-2;
  c.velocity_iterations = 12;
  c.pressure_tolerance = 1e-3;
  c.pressure_iterations = 60;
  c.eta0 = 1e21;
  c.velocity_scale = 1;
  c.density = 3300;
  c.anomaly = 80;
  c.gravity[1] = 9.80665;
  c.width = .12;
  c.contrast = 10;
  c.eta_min = 1e18;
  c.eta_max = 1e25;
  c.viscosity = 1; /* 0 constant, 1 linear, 2 inclusion. Mode 0 manufactured. */
  for (f = 0; f < 6; f++)
    c.bc[f].kind = EXACT;
  return c;
}
static double number(const char *text, int allow_nan) {
  char *end;
  double v = strtod(text, &end);
  if (end == text || *end || (!isfinite(v) && !(allow_nan && isnan(v))))
    fail("Invalid number in configuration");
  return v;
}
static int integer(const char *text) {
  double v = number(text, 0);
  if (v != floor(v) || v < 0 || v > 100000)
    fail("Expected a nonnegative integer");
  return (int)v;
}
static void read_config(Config *c, const char *path) {
  FILE *f = fopen(path, "r");
  char line[1024];
  if (!f)
    fail("Cannot open input configuration");
  while (fgets(line, sizeof(line), f)) {
    char *t[8], *p;
    int n = 0, d, face = -1;
    p = strchr(line, '#');
    if (p)
      *p = 0;
    p = strtok(line, " \t\r\n=");
    while (p && n < 8) {
      t[n++] = p;
      p = strtok(NULL, " \t\r\n=");
    }
    if (n == 0)
      continue;
    if (!strcmp(t[0], "lengths_km") || !strcmp(t[0], "elements") ||
        !strcmp(t[0], "velocity_nodes") || !strcmp(t[0], "gravity_m_s2") ||
        !strcmp(t[0], "inclusion_center")) {
      if (n != 4)
        fail("Vector input needs three entries");
      for (d = 0; d < 3; d++) {
        if (!strcmp(t[0], "lengths_km"))
          c->length[d] = number(t[d + 1], 0);
        else if (!strcmp(t[0], "elements"))
          c->n[d] = integer(t[d + 1]);
        else if (!strcmp(t[0], "velocity_nodes"))
          c->requested_nodes[d] = integer(t[d + 1]);
        else if (!strcmp(t[0], "gravity_m_s2"))
          c->gravity[d] = number(t[d + 1], 0);
        else
          c->center[d] = number(t[d + 1], 0);
      }
      continue;
    }
    for (d = 0; d < 6; d++)
      if (!strcmp(t[0], face_names[d]))
        face = d;
    if (face >= 0) {
      Boundary *b = &c->bc[face];
      if (n < 2)
        fail("Missing boundary type");
      if (!strcmp(t[1], "exact"))
        b->kind = EXACT;
      else if (!strcmp(t[1], "noslip"))
        b->kind = NOSLIP;
      else if (!strcmp(t[1], "free_slip"))
        b->kind = FREESLIP;
      else if (!strcmp(t[1], "velocity"))
        b->kind = VELOCITY;
      else if (!strcmp(t[1], "traction"))
        b->kind = TRACTION;
      else
        fail("Unknown boundary type");
      if (b->kind >= VELOCITY) {
        if (n != 5)
          fail("velocity/traction needs three components");
        for (d = 0; d < 3; d++)
          b->value[d] = number(t[d + 2], b->kind == VELOCITY);
      } else if (n != 2)
        fail("This boundary type takes no additional numbers");
      continue;
    }
    if (n != 2)
      fail("Scalar configuration entry needs one value");
    if (!strcmp(t[0], "mesh_input")) {
      if (!strcmp(t[1], "elements"))
        c->mesh_nodes = 0;
      else if (!strcmp(t[1], "nodes"))
        c->mesh_nodes = 1;
      else
        fail("mesh_input must be elements or nodes");
    } else if (!strcmp(t[0], "relative_tolerance"))
      c->tolerance = number(t[1], 0);
    else if (!strcmp(t[0], "restart"))
      c->restart = integer(t[1]);
    else if (!strcmp(t[0], "max_iterations"))
      c->max_iterations = integer(t[1]);
    else if (!strcmp(t[0], "velocity_inner_tolerance"))
      c->velocity_tolerance = number(t[1], 0);
    else if (!strcmp(t[0], "velocity_inner_iterations"))
      c->velocity_iterations = integer(t[1]);
    else if (!strcmp(t[0], "pressure_inner_tolerance"))
      c->pressure_tolerance = number(t[1], 0);
    else if (!strcmp(t[0], "pressure_inner_iterations"))
      c->pressure_iterations = integer(t[1]);
    else if (!strcmp(t[0], "velocity_degree"))
      c->degree = integer(t[1]);
    else if (!strcmp(t[0], "benchmark")) {
      if (!strcmp(t[1], "polynomial"))
        c->benchmark = 0;
      else if (!strcmp(t[1], "trigonometric"))
        c->benchmark = 1;
      else
        fail("Unknown benchmark");
    } else if (!strcmp(t[0], "mode")) {
      if (!strcmp(t[1], "manufactured"))
        c->mode = 0;
      else if (!strcmp(t[1], "gravity"))
        c->mode = 1;
      else if (!strcmp(t[1], "custom"))
        c->mode = 2;
      else
        fail("Unknown mode");
    } else if (!strcmp(t[0], "viscosity_model")) {
      if (!strcmp(t[1], "constant"))
        c->viscosity = 0;
      else if (!strcmp(t[1], "linear"))
        c->viscosity = 1;
      else if (!strcmp(t[1], "inclusion"))
        c->viscosity = 2;
      else
        fail("Unknown viscosity model");
    } else if (!strcmp(t[0], "averaging")) {
      if (!strcmp(t[1], "arithmetic"))
        c->averaging = 0;
      else if (!strcmp(t[1], "geometric"))
        c->averaging = 1;
      else if (!strcmp(t[1], "harmonic"))
        c->averaging = 2;
      else
        fail("Unknown averaging");
    } else if (!strcmp(t[0], "quadrature"))
      c->order = integer(t[1]);
    else if (!strcmp(t[0], "error_quadrature"))
      c->error_order = integer(t[1]);
    else if (!strcmp(t[0], "max_equations"))
      c->max_equations = integer(t[1]);
    else if (!strcmp(t[0], "eta0_Pa_s"))
      c->eta0 = number(t[1], 0);
    else if (!strcmp(t[0], "velocity_scale_cm_year"))
      c->velocity_scale = number(t[1], 0);
    else if (!strcmp(t[0], "pressure_mean_Pa"))
      c->pressure_mean = number(t[1], 0);
    else if (!strcmp(t[0], "grad_div"))
      c->gamma = number(t[1], 0);
    else if (!strcmp(t[0], "density_kg_m3"))
      c->density = number(t[1], 0);
    else if (!strcmp(t[0], "density_anomaly_kg_m3"))
      c->anomaly = number(t[1], 0);
    else if (!strcmp(t[0], "inclusion_width"))
      c->width = number(t[1], 0);
    else if (!strcmp(t[0], "viscosity_contrast"))
      c->contrast = number(t[1], 0);
    else if (!strcmp(t[0], "eta_min_Pa_s"))
      c->eta_min = number(t[1], 0);
    else if (!strcmp(t[0], "eta_max_Pa_s"))
      c->eta_max = number(t[1], 0);
    else
      fail("Unknown configuration key");
  }
  fclose(f);
  if (c->degree < 2 || c->degree > 4)
    fail("velocity_degree must be 2, 3, or 4");
  if (c->mesh_nodes)
    for (int d = 0; d < 3; d++) {
      int nodes = c->requested_nodes[d];
      if (nodes < 2 * c->degree + 1 || (nodes - 1) % c->degree != 0)
        fail("For each direction, velocity_nodes must equal degree*elements+1, "
             "with elements >= 2");
      c->n[d] = (nodes - 1) / c->degree;
    }
  if (c->tolerance <= 0 || c->tolerance >= 1 || c->restart < 2 ||
      c->restart > 300 || c->max_iterations < 1)
    fail("Require 0 < relative_tolerance < 1, restart 2..300, max_iterations "
         "1..100000");
  if (c->velocity_tolerance <= 0 || c->velocity_tolerance >= 1 ||
      c->pressure_tolerance <= 0 || c->pressure_tolerance >= 1 ||
      c->velocity_iterations < 1 || c->pressure_iterations < 1)
    fail("Inner tolerances must be between zero and one; inner iteration "
         "limits must be positive");
  if (c->order < c->degree + 1 || c->error_order < c->degree + 2)
    fail("Use quadrature >= degree+1 and error_quadrature >= degree+2");
  for (int d = 0; d < 3; d++)
    if (c->n[d] < 2 || c->n[d] > 100 || c->length[d] <= 0 || c->center[d] < 0 ||
        c->center[d] > 1)
      fail("Invalid lengths, element counts or inclusion center");
  if (c->order < 3 || c->order > QMAX || c->error_order < 3 ||
      c->error_order > QMAX)
    fail("Quadrature must be 3 to 10");
  if (c->eta0 <= 0 || c->velocity_scale <= 0 || c->width <= 0 ||
      c->contrast <= 0 || c->density <= 0 ||
      c->density + fmin(0, c->anomaly) <= 0 || c->gamma < 0)
    fail("Invalid positive physical parameter");
  if (c->eta_min <= 0 || c->eta_max < c->eta_min || c->max_equations < 100 ||
      c->max_equations > 100000)
    fail("Invalid viscosity bounds or equation cap (100..100000)");
  if (c->mode == 0 && (c->viscosity == 2 || c->eta_min > c->eta0 ||
                       c->eta_max < c->eta0 * (c->viscosity == 1 ? 4 : 1)))
    fail("Manufactured test requires unclipped constant/linear viscosity");
  for (int f0 = 0; f0 < 6; f0++)
    if (c->mode != 0 && c->bc[f0].kind == EXACT)
      fail("Exact walls require manufactured mode");
}
/* Tensor-product Gauss rule: Newton iteration on Legendre polynomials. */
static void gauss(int n, double *x, double *w) {
  int i, k, j;
  for (i = 0; i < (n + 1) / 2; i++) {
    double z = cos(3.14159265358979323846 * (i + .75) / (n + .5)), pp = 0;
    for (k = 0; k < 50; k++) {
      double p1 = 1, p2 = 0, p3, old;
      for (j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
      }
      pp = n * (z * p1 - p2) / (z * z - 1);
      old = z;
      z -= p1 / pp;
      if (fabs(z - old) < 2e-15)
        break;
    }
    if (k == 50)
      fail("Gauss quadrature root failed to converge");
    x[i] = -z;
    x[n - 1 - i] = z;
    w[i] = w[n - 1 - i] = 2 / ((1 - z * z) * pp * pp);
  }
}
/* Local order is x then y then z, with z changing fastest. */

/* Lagrange interpolation on equally spaced reference points in [-1,1].
 * The derivative is computed by differentiating the product explicitly.
 */
static void lagrange(int degree, double t, double *value, double *derivative) {
  for (int a = 0; a <= degree; a++) {
    double ta = -1 + 2.0 * a / degree;
    value[a] = 1;
    derivative[a] = 0;
    for (int b = 0; b <= degree; b++)
      if (b != a) {
        double tb = -1 + 2.0 * b / degree;
        value[a] *= (t - tb) / (ta - tb);
      }
    for (int k = 0; k <= degree; k++)
      if (k != a) {
        double tk = -1 + 2.0 * k / degree, term = 1 / (ta - tk);
        for (int b = 0; b <= degree; b++)
          if (b != a && b != k) {
            double tb = -1 + 2.0 * b / degree;
            term *= (t - tb) / (ta - tb);
          }
        derivative[a] += term;
      }
  }
}
static void basis(const Model *model, const double xi[3], const double h[3],
                  double N[MAXV], double D[MAXV][3], double M[MAXP]) {
  int degree = model->c.degree, a = 0;
  double l[3][5], dl[3][5], p[3][5], unused[5];
  for (int d = 0; d < 3; d++) {
    lagrange(degree, xi[d], l[d], dl[d]);
    lagrange(degree - 1, xi[d], p[d], unused);
    for (int i = 0; i <= degree; i++)
      dl[d][i] *= 2 / h[d];
  }
  for (int i = 0; i <= degree; i++)
    for (int j = 0; j <= degree; j++)
      for (int k = 0; k <= degree; k++, a++) {
        N[a] = l[0][i] * l[1][j] * l[2][k];
        D[a][0] = dl[0][i] * l[1][j] * l[2][k];
        D[a][1] = l[0][i] * dl[1][j] * l[2][k];
        D[a][2] = l[0][i] * l[1][j] * dl[2][k];
      }
  a = 0;
  for (int i = 0; i < degree; i++)
    for (int j = 0; j < degree; j++)
      for (int k = 0; k < degree; k++, a++)
        M[a] = p[0][i] * p[1][j] * p[2][k];
}
static void connectivity(const Model *m, const int cell[3], int v[MAXV],
                         int p[MAXP]) {
  int degree = m->c.degree, a = 0;
  for (int i = 0; i <= degree; i++)
    for (int j = 0; j <= degree; j++)
      for (int k = 0; k <= degree; k++, a++)
        v[a] = node_id(degree * cell[0] + i, degree * cell[1] + j,
                       degree * cell[2] + k, m->vn);
  a = 0;
  for (int i = 0; i < degree; i++)
    for (int j = 0; j < degree; j++)
      for (int k = 0; k < degree; k++, a++)
        p[a] = node_id((degree - 1) * cell[0] + i, (degree - 1) * cell[1] + j,
                       (degree - 1) * cell[2] + k, m->pn);
}
/* Derivatives of the same curl-based manufactured field used by the notebook.
 */
static double derivative(const Model *m, double poly[4][3], int component,
                         const int extra[3]) {
  static const double coeff[3][2] = {{3, -2}, {1, -3}, {2, -1}};
  static const int orders[3][2][3] = {
      {{0, 1, 0}, {0, 0, 1}}, {{0, 0, 1}, {1, 0, 0}}, {{1, 0, 0}, {0, 1, 0}}};
  double sum = 0;
  int t, d;
  for (t = 0; t < 2; t++) {
    double v = 256 * coeff[component][t] * m->ell[component];
    for (d = 0; d < 3; d++)
      v *= poly[orders[component][t][d] + extra[d]][d] /
           pow(m->ell[d], extra[d]);
    sum += v;
  }
  return sum;
}
static void exact(const Model *m, const double x[3], double u[3],
                  double G[3][3], double lap[3], double *p, double gp[3]) {
  double poly[4][3];
  int i, j, extra[3] = {0, 0, 0};
  *p = -.75;
  for (i = 0; i < 3; i++) {
    double t = x[i] / m->ell[i];
    poly[0][i] = t * t * (1 - t) * (1 - t);
    poly[1][i] = 2 * t - 6 * t * t + 4 * t * t * t;
    poly[2][i] = 2 - 12 * t + 12 * t * t;
    poly[3][i] = -12 + 24 * t;
    *p += t * t * t;
    gp[i] = 3 * t * t / m->ell[i];
    if (m->c.benchmark == 1) {
      const double pi = 3.14159265358979323846;
      poly[0][i] = sin(pi * t) * sin(pi * t) / 16;
      poly[1][i] = pi * sin(2 * pi * t) / 16;
      poly[2][i] = 2 * pi * pi * cos(2 * pi * t) / 16;
      poly[3][i] = -4 * pi * pi * pi * sin(2 * pi * t) / 16;
    }
  }
  for (i = 0; i < 3; i++) {
    u[i] = derivative(m, poly, i, extra);
    lap[i] = 0;
    for (j = 0; j < 3; j++) {
      extra[j] = 1;
      G[i][j] = derivative(m, poly, i, extra);
      extra[j] = 2;
      lap[i] += derivative(m, poly, i, extra);
      extra[j] = 0;
    }
  }
}
/* Edit this C function for a custom physical force. Input metres; output N/m^3.
 */
static void custom_force(const double x_m[3], double force[3]) {
  (void)x_m;
  force[0] = 0;
  force[1] = 0;
  force[2] = 0;
}
static double coefficients(const Model *m, const int cell[3], const double x[3],
                           const double M[MAXP], double f[3]) {
  const Config *c = &m->c;
  (void)M;
  double eta = 1, rho = 0, acc = 0;
  int i, j, k, a = 0, d;
  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++, a++) {
        double r[3] = {(cell[0] + i) / (double)c->n[0],
                       (cell[1] + j) / (double)c->n[1],
                       (cell[2] + k) / (double)c->n[2]},
               dist = 0, blob, v;
        for (d = 0; d < 3; d++)
          dist += (r[d] - c->center[d]) * (r[d] - c->center[d]);
        blob = exp(-dist / (c->width * c->width));
        v = exp(log(c->contrast) * blob);
        double local[3];
        for (int d0 = 0; d0 < 3; d0++)
          local[d0] = x[d0] / m->h[d0] - cell[d0];
        double material_weight = (i ? local[0] : 1 - local[0]) *
                                 (j ? local[1] : 1 - local[1]) *
                                 (k ? local[2] : 1 - local[2]);
        rho += material_weight * (c->density + c->anomaly * blob);
        acc += material_weight * (c->averaging == 0   ? v
                                  : c->averaging == 1 ? log(v)
                                                      : 1 / v);
      }
  if (c->viscosity == 1)
    for (d = 0; d < 3; d++)
      eta += x[d] / m->ell[d];
  if (c->viscosity == 2)
    eta = c->averaging == 0 ? acc : c->averaging == 1 ? exp(acc) : 1 / acc;
  eta = fmax(c->eta_min / c->eta0, fmin(c->eta_max / c->eta0, eta));
  if (c->mode == 0) {
    double u[3], G[3][3], lap[3], p, gp[3];
    exact(m, x, u, G, lap, &p, gp);
    for (i = 0; i < 3; i++) {
      f[i] = -eta * lap[i] + gp[i];
      if (c->viscosity == 1)
        for (j = 0; j < 3; j++)
          f[i] -= (G[i][j] + G[j][i]) / m->ell[j];
    }
  } else if (c->mode == 1)
    for (i = 0; i < 3; i++)
      f[i] = rho * c->gravity[i] / m->F0;
  else {
    double xm[3];
    for (i = 0; i < 3; i++)
      xm[i] = x[i] * m->L0;
    custom_force(xm, f);
    for (i = 0; i < 3; i++)
      f[i] /= m->F0;
  }
  for (i = 0; i < 3; i++)
    if (!isfinite(f[i]) || !isfinite(eta))
      fail("Nonfinite material or force");
  return eta;
}
/* Apply essential velocity conditions. NaN leaves a velocity component natural.
 */
void boundaries(Model *m) {
  int i, j, k, d, face, id, free_count = 0;
  double rigid_gram[6][6] = {{0}};
  for (id = 0; id < 3 * m->nv; id++)
    m->map[id] = -2;
  for (i = 0; i < m->vn[0]; i++)
    for (j = 0; j < m->vn[1]; j++)
      for (k = 0; k < m->vn[2]; k++) {
        int index[3] = {i, j, k};
        double x[3] = {i * m->h[0] / m->c.degree, j * m->h[1] / m->c.degree,
                       k * m->h[2] / m->c.degree};
        id = node_id(i, j, k, m->vn);
        for (face = 0; face < 6; face++) {
          Boundary b = m->c.bc[face];
          int axis = face / 2, side = face % 2;
          if (index[axis] != (side ? m->vn[axis] - 1 : 0))
            continue;
          for (d = 0; d < 3; d++) {
            int fixed = b.kind == EXACT || b.kind == NOSLIP ||
                        (b.kind == FREESLIP && d == axis) ||
                        (b.kind == VELOCITY && !isnan(b.value[d]));
            double value =
                b.kind == VELOCITY ? b.value[d] / m->c.velocity_scale : 0;
            if (!fixed)
              continue;
            if (m->map[3 * id + d] == -1 &&
                fabs(m->ub[3 * id + d] - value) > 1e-10)
              fail("Conflicting prescribed velocities at an edge or corner");
            m->map[3 * id + d] = -1;
            m->ub[3 * id + d] = value;
          }
        }
        for (d = 0; d < 3; d++)
          if (m->map[3 * id + d] == -1) {
            double row[6] = {0};
            row[d] = 1;
            if (d == 0) {
              row[4] = x[2];
              row[5] = -x[1];
            }
            if (d == 1) {
              row[3] = -x[2];
              row[5] = x[0];
            }
            if (d == 2) {
              row[3] = x[1];
              row[4] = -x[0];
            }
            for (int a = 0; a < 6; a++)
              for (int b = 0; b < 6; b++)
                rigid_gram[a][b] += row[a] * row[b];
          }
      }
  /* Rank check on constrained rigid-body modes prevents a singular physical
   * model. */
  for (i = 0; i < 6; i++) {
    int pivot = i;
    for (j = i + 1; j < 6; j++)
      if (fabs(rigid_gram[j][i]) > fabs(rigid_gram[pivot][i]))
        pivot = j;
    if (fabs(rigid_gram[pivot][i]) < 1e-12)
      fail("Boundaries leave an unconstrained rigid motion");
    for (k = i; k < 6; k++) {
      double t = rigid_gram[i][k];
      rigid_gram[i][k] = rigid_gram[pivot][k];
      rigid_gram[pivot][k] = t;
    }
    for (j = i + 1; j < 6; j++) {
      double t = rigid_gram[j][i] / rigid_gram[i][i];
      for (k = i; k < 6; k++)
        rigid_gram[j][k] -= t * rigid_gram[i][k];
    }
  }
  for (id = 0; id < 3 * m->nv; id++)
    if (m->map[id] != -1)
      m->map[id] = free_count++;
  m->nf = free_count;
  if (!free_count)
    fail("No unknown velocity coefficients");
}
void traction(const Model *m, const int cell[3], const double *q,
              const double *w, double F[MAXDOF]) {
  int face, a, b, n, i;
  for (face = 0; face < 6; face++) {
    int axis = face / 2, side = face % 2, other[2], s = 0;
    if (m->c.bc[face].kind != TRACTION ||
        cell[axis] != (side ? m->c.n[axis] - 1 : 0))
      continue;
    for (i = 0; i < 3; i++)
      if (i != axis)
        other[s++] = i;
    for (a = 0; a < m->c.order; a++)
      for (b = 0; b < m->c.order; b++) {
        double xi[3], N[MAXV], D[MAXV][3], M[MAXP],
            W = w[a] * w[b] * m->h[other[0]] * m->h[other[1]] / 4;
        xi[axis] = side ? 1 : -1;
        xi[other[0]] = q[a];
        xi[other[1]] = q[b];
        basis(m, xi, m->h, N, D, M);
        for (n = 0; n < m->local_v; n++)
          for (i = 0; i < 3; i++)
            F[3 * n + i] += W * N[n] * m->c.bc[face].value[i] / m->P0;
      }
  }
}
static void evaluate(const Model *m, const double x[3], double u[3],
                     double G[3][3], double *p) {
  int cell[3], v[MAXV], pn[MAXP], d, a;
  double xi[3], N[MAXV], D[MAXV][3], M[MAXP];
  for (d = 0; d < 3; d++) {
    cell[d] = (int)floor(x[d] / m->h[d]);
    if (cell[d] < 0)
      cell[d] = 0;
    if (cell[d] >= m->c.n[d])
      cell[d] = m->c.n[d] - 1;
    xi[d] = 2 * (x[d] / m->h[d] - cell[d]) - 1;
  }
  basis(m, xi, m->h, N, D, M);
  connectivity(m, cell, v, pn);
  memset(u, 0, 3 * sizeof(double));
  memset(G, 0, 9 * sizeof(double));
  *p = 0;
  for (a = 0; a < m->local_v; a++)
    for (d = 0; d < 3; d++) {
      u[d] += N[a] * m->velocity[3 * v[a] + d];
      for (int j = 0; j < 3; j++)
        G[d][j] += D[a][j] * m->velocity[3 * v[a] + d];
    }
  for (a = 0; a < m->local_p; a++)
    *p += M[a] * m->pressure[pn[a]];
}
static int has_exact(const Model *m) {
  if (m->c.mode != 0)
    return 0;
  for (int f = 0; f < 6; f++)
    if (m->c.bc[f].kind != EXACT && m->c.bc[f].kind != NOSLIP)
      return 0;
  return 1;
}
void errors_and_stats(const Model *m, const char *prefix) {
  double q[QMAX], w[QMAX], sum[14] = {0}, maximum[3] = {0}, divergence = 0;
  int cell[3], o = m->c.error_order, known = has_exact(m);
  const char *names[14] = {"velocity_L2",
                           "pressure_L2",
                           "velocity_H1_semi",
                           "speed_L2",
                           "strain_L2",
                           "viscous_stress_L2",
                           "total_stress_L2",
                           "energy_error",
                           "exact_velocity_norm",
                           "exact_pressure_norm",
                           "vx_L2",
                           "vy_L2",
                           "vz_L2",
                           "unused"};
  gauss(o, q, w);
  for (cell[0] = 0; cell[0] < m->c.n[0]; cell[0]++)
    for (cell[1] = 0; cell[1] < m->c.n[1]; cell[1]++)
      for (cell[2] = 0; cell[2] < m->c.n[2]; cell[2]++)
        for (int ix = 0; ix < o; ix++)
          for (int iy = 0; iy < o; iy++)
            for (int iz = 0; iz < o; iz++) {
              double xi[3] = {q[ix], q[iy], q[iz]}, x[3], u[3], G[3][3], p,
                     ue[3], Ge[3][3], pe, lap[3], gp[3];
              double W =
                  w[ix] * w[iy] * w[iz] * m->h[0] * m->h[1] * m->h[2] / 8;
              for (int d = 0; d < 3; d++)
                x[d] = (cell[d] + (xi[d] + 1) / 2) * m->h[d];
              evaluate(m, x, u, G, &p);
              double div = G[0][0] + G[1][1] + G[2][2];
              divergence += W * div * div;
              if (known) {
                double N[MAXV], D[MAXV][3], M[MAXP], force[3], eta,
                    ev = 0, eg = 0, strain = 0, stress = 0, totalstress = 0, dp,
                    speed;
                basis(m, xi, m->h, N, D, M);
                eta = coefficients(m, cell, x, M, force);
                exact(m, x, ue, Ge, lap, &pe, gp);
                pe += m->c.pressure_mean / m->P0;
                dp = p - pe;
                speed = norm(u, 3) - norm(ue, 3);
                for (int i = 0; i < 3; i++) {
                  double du = u[i] - ue[i];
                  ev += du * du;
                  sum[10 + i] += W * du * du;
                  for (int j = 0; j < 3; j++) {
                    double dg = G[i][j] - Ge[i][j],
                           e = .5 * (dg + G[j][i] - Ge[j][i]), s = 2 * eta * e,
                           t = s - (i == j ? dp : 0);
                    eg += dg * dg;
                    strain += e * e;
                    stress += s * s;
                    totalstress += t * t;
                  }
                }
                double values[10] = {ev,
                                     dp * dp,
                                     eg,
                                     speed * speed,
                                     strain,
                                     stress,
                                     totalstress,
                                     2 * eta * strain,
                                     norm(ue, 3) * norm(ue, 3),
                                     pe * pe};
                for (int i = 0; i < 10; i++)
                  sum[i] += W * values[i];
                maximum[0] = fmax(maximum[0], sqrt(ev));
                maximum[1] = fmax(maximum[1], fabs(dp));
                maximum[2] = fmax(maximum[2], fabs(speed));
              }
            }
  FILE *f = output_file(prefix, "_stats.json");
  double pressure_mean = 0, volume = 0;
  for (int i = 0; i < m->np; i++) {
    pressure_mean += m->mean_weights[i] * m->pressure[i];
    volume += m->mean_weights[i];
  }
  fprintf(f,
          "{\n\"velocity_nodes\":%d,\n\"pressure_nodes\":%d,\n\"elements\":%d,"
          "\n\"equations\":%d,"
          "\n\"pressure_gauge\":%s,\n\"has_exact\":%s,\n",
          m->nv, m->np, m->c.n[0] * m->c.n[1] * m->c.n[2], m->neq,
          m->gauge ? "true" : "false", known ? "true" : "false");
  fprintf(f,
          "\"scaled_residual\":%.17g,\n\"continuity_residual\":%.17g,"
          "\n\"boundary_error\":0,"
          "\n\"pressure_mean_Pa\":%.17g,\n\"gauge_multiplier\":%.17g,"
          "\n\"matrix_symmetry_error\":"
          "%.17g,\n",
          m->residual, m->continuity, pressure_mean / volume * m->P0, m->lambda,
          m->symmetry);
  fprintf(f, "\"velocity_degree\":%d,\n\"pressure_degree\":%d,\n", m->c.degree,
          m->c.degree - 1);
  fprintf(f, "\"benchmark\":\"%s\",\n\"velocity_scale_cm_year\":%.17g,\n",
          m->c.benchmark ? "trigonometric" : "polynomial", m->c.velocity_scale);
  fprintf(f,
          "\"assembly_seconds\":%.17g,\n\"solve_seconds\":%.17g,\n\"sparse_"
          "matrix_storage_"
          "MiB\":%.17g,\n",
          m->assembly_time, m->solve_time, m->sparse_storage_MiB);
  if (known) {
    for (int i = 0; i < 13; i++)
      fprintf(f, "\"%s\":%.17g,\n", names[i], sqrt(sum[i]));
    fprintf(f,
            "\"velocity_relative_L2\":%.17g,\n\"pressure_relative_L2\":%.17g,"
            "\n\"velocity_Linf_"
            "sampled\":%.17g,\n\"pressure_Linf_sampled\":%.17g,\n\"speed_Linf_"
            "sampled\":%.17g,\n",
            sqrt(sum[0] / sum[8]), sqrt(sum[1] / sum[9]), maximum[0],
            maximum[1], maximum[2]);
  }
  fprintf(f, "\"divergence_L2\":%.17g\n}\n", sqrt(divergence));
  fclose(f);
  printf("nodes v=%d p=%d; equations=%d; residual=%.3e; assembly=%.3fs "
         "solve=%.3fs\n",
         m->nv, m->np, m->neq, m->residual, m->assembly_time, m->solve_time);
  if (known)
    printf("relative L2 errors: velocity %.10g%%; pressure %.10g%%\n",
           100 * sqrt(sum[0] / sum[8]), 100 * sqrt(sum[1] / sum[9]));
}

/* Physical-unit field exports. Derivatives are evaluated from FEM shape
   functions, never from finite differences of plotted samples. Stress may jump
   across cells. */
static void all_fields(const Model *m, const double x[3], double v[17],
                       double e[17]) {
  double u[3], G[3][3], p, ue[3], Ge[3][3], pe, lap[3], gp[3];
  double xi[3], N[MAXV], D[MAXV][3], M[MAXP], force[3], eta;
  int cell[3];
  const int a[6] = {0, 1, 2, 0, 0, 1}, b[6] = {0, 1, 2, 1, 2, 2};
  for (int d = 0; d < 3; d++) {
    cell[d] = (int)floor(x[d] / m->h[d]);
    if (cell[d] >= m->c.n[d])
      cell[d] = m->c.n[d] - 1;
    xi[d] = 2 * (x[d] / m->h[d] - cell[d]) - 1;
  }
  basis(m, xi, m->h, N, D, M);
  eta = coefficients(m, cell, x, M, force);
  evaluate(m, x, u, G, &p);
  exact(m, x, ue, Ge, lap, &pe, gp);
  pe += m->c.pressure_mean / m->P0;
  for (int d = 0; d < 3; d++) {
    v[d] = u[d] * m->c.velocity_scale;
    e[d] = ue[d] * m->c.velocity_scale;
  }
  v[3] = norm(u, 3) * m->c.velocity_scale;
  e[3] = norm(ue, 3) * m->c.velocity_scale;
  v[4] = p * m->P0 / 1e6;
  e[4] = pe * m->P0 / 1e6;
  for (int k = 0; k < 6; k++) {
    int i = a[k], j = b[k];
    v[5 + k] = eta * (G[i][j] + G[j][i]) * m->P0 / 1e6;
    e[5 + k] = eta * (Ge[i][j] + Ge[j][i]) * m->P0 / 1e6;
    v[11 + k] = v[5 + k] - (i == j ? v[4] : 0);
    e[11 + k] = e[5 + k] - (i == j ? e[4] : 0);
  }
  if (!has_exact(m))
    for (int k = 0; k < 17; k++)
      e[k] = NAN;
}
static void field_header(FILE *f) {
  const char *names[17] = {
      "vx",         "vy",         "vz",         "speed",      "pressure",
      "viscous_xx", "viscous_yy", "viscous_zz", "viscous_xy", "viscous_xz",
      "viscous_yz", "total_xx",   "total_yy",   "total_zz",   "total_xy",
      "total_xz",   "total_yz"};
  fprintf(f, "x_km,y_km,z_km");
  for (int k = 0; k < 17; k++)
    fprintf(f, ",%s,exact_%s", names[k], names[k]);
  fprintf(f, "\n");
}
static void field_row(FILE *f, const Model *m, const double x[3]) {
  double v[17], e[17];
  all_fields(m, x, v, e);
  fprintf(f, "%.17g,%.17g,%.17g", x[0] * m->L0 / 1000, x[1] * m->L0 / 1000,
          x[2] * m->L0 / 1000);
  for (int k = 0; k < 17; k++)
    fprintf(f, ",%.17g,%.17g", v[k], e[k]);
  fprintf(f, "\n");
}
void export_all_fields(const Model *m, const char *prefix) {
  const char *suffix[3] = {"_fields_x.csv", "_fields_y.csv", "_fields_z.csv"};
  for (int axis = 0; axis < 3; axis++) {
    FILE *f = output_file(prefix, suffix[axis]);
    field_header(f);
    for (int i = 0; i <= 200; i++) {
      double x[3] = {.31 * m->ell[0], .37 * m->ell[1], .43 * m->ell[2]};
      x[axis] = i * m->ell[axis] / 200;
      field_row(f, m, x);
    }
    fclose(f);
  }
  FILE *f = output_file(prefix, "_fields_slice.csv");
  field_header(f);
  for (int i = 0; i <= 40; i++)
    for (int j = 0; j <= 40; j++) {
      double x[3] = {i * m->ell[0] / 40, j * m->ell[1] / 40, .43 * m->ell[2]};
      field_row(f, m, x);
    }
  fclose(f);
  f = output_file(prefix, "_fields_volume.csv");
  field_header(f);
  for (int i = 0; i <= 16; i++)
    for (int j = 0; j <= 16; j++)
      for (int k = 0; k <= 16; k++) {
        double x[3] = {i * m->ell[0] / 16, j * m->ell[1] / 16,
                       k * m->ell[2] / 16};
        field_row(f, m, x);
      }
  fclose(f);
}

void export_results(const Model *m, const char *prefix) {
  FILE *f = output_file(prefix, "_velocity.csv");
  int i, j, k, d;
  fprintf(
      f, "x_km,y_km,z_km,vx_dimensionless,vy_dimensionless,vz_dimensionless\n");
  for (i = 0; i < m->vn[0]; i++)
    for (j = 0; j < m->vn[1]; j++)
      for (k = 0; k < m->vn[2]; k++) {
        int id = node_id(i, j, k, m->vn);
        fprintf(f, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                i * m->c.length[0] / (m->vn[0] - 1),
                j * m->c.length[1] / (m->vn[1] - 1),
                k * m->c.length[2] / (m->vn[2] - 1), m->velocity[3 * id],
                m->velocity[3 * id + 1], m->velocity[3 * id + 2]);
      }
  fclose(f);
  f = output_file(prefix, "_pressure.csv");
  fprintf(f, "x_km,y_km,z_km,pressure_dimensionless\n");
  for (i = 0; i < m->pn[0]; i++)
    for (j = 0; j < m->pn[1]; j++)
      for (k = 0; k < m->pn[2]; k++)
        fprintf(f, "%.17g,%.17g,%.17g,%.17g\n",
                i * m->c.length[0] / (m->pn[0] - 1),
                j * m->c.length[1] / (m->pn[1] - 1),
                k * m->c.length[2] / (m->pn[2] - 1),
                m->pressure[node_id(i, j, k, m->pn)]);
  fclose(f);
  f = output_file(prefix, "_slice.csv");
  fprintf(f, "x_km,y_km,z_km,vx_cm_year,vy_cm_year,vz_cm_year,speed_cm_year,"
             "pressure_MPa,exact_vx_"
             "cm_year,exact_vy_cm_year,exact_vz_cm_year,exact_speed_cm_year,"
             "exact_pressure_MPa\n");
  for (i = 0; i < 41; i++)
    for (j = 0; j < 41; j++) {
      double x[3] = {i * m->ell[0] / 40, j * m->ell[1] / 40, .43 * m->ell[2]},
             u[3], G[3][3], p, ue[3], Ge[3][3], lap[3], pe, gp[3];
      evaluate(m, x, u, G, &p);
      exact(m, x, ue, Ge, lap, &pe, gp);
      pe += m->c.pressure_mean / m->P0;
      fprintf(f, "%.17g,%.17g,%.17g", x[0] * m->L0 / 1000, x[1] * m->L0 / 1000,
              x[2] * m->L0 / 1000);
      for (d = 0; d < 3; d++)
        fprintf(f, ",%.17g", u[d] * m->c.velocity_scale);
      fprintf(f, ",%.17g,%.17g", norm(u, 3) * m->c.velocity_scale,
              p * m->P0 / 1e6);
      if (has_exact(m)) {
        for (d = 0; d < 3; d++)
          fprintf(f, ",%.17g", ue[d] * m->c.velocity_scale);
        fprintf(f, ",%.17g,%.17g\n", norm(ue, 3) * m->c.velocity_scale,
                pe * m->P0 / 1e6);
      } else
        fprintf(f, ",nan,nan,nan,nan,nan\n");
    }
  fclose(f);

  f = output_file(prefix, "_profile.csv");
  fprintf(f, "x_km,vx,vy,vz,speed,pressure,exact_vx,exact_vy,exact_vz,exact_"
             "speed,exact_pressure\n");
  for (i = 0; i < 201; i++) {
    double x[3] = {i * m->ell[0] / 200, .37 * m->ell[1], .43 * m->ell[2]}, u[3],
           G[3][3], p, ue[3], Ge[3][3], lap[3], pe, gp[3];
    evaluate(m, x, u, G, &p);
    exact(m, x, ue, Ge, lap, &pe, gp);
    pe += m->c.pressure_mean / m->P0;
    fprintf(f, "%.17g", x[0] * m->L0 / 1000);
    for (d = 0; d < 3; d++)
      fprintf(f, ",%.17g", u[d] * m->c.velocity_scale);
    fprintf(f, ",%.17g,%.17g", norm(u, 3) * m->c.velocity_scale,
            p * m->P0 / 1e6);
    if (has_exact(m)) {
      for (d = 0; d < 3; d++)
        fprintf(f, ",%.17g", ue[d] * m->c.velocity_scale);
      fprintf(f, ",%.17g,%.17g\n", norm(ue, 3) * m->c.velocity_scale,
              pe * m->P0 / 1e6);
    } else
      fprintf(f, ",nan,nan,nan,nan,nan\n");
  }
  fclose(f);
  /* Legacy VTK uses x-fastest ordering; reorder explicitly for 3-D viewers. */
  f = output_file(prefix, ".vtk");
  fprintf(f,
          "# vtk DataFile Version 3.0\nQp/Q(p-1) Stokes; coordinates km, "
          "velocity cm/year, pressure "
          "MPa\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS %d %d %d\nPOINTS %d "
          "double\n",
          m->vn[0], m->vn[1], m->vn[2], m->nv);
  for (k = 0; k < m->vn[2]; k++)
    for (j = 0; j < m->vn[1]; j++)
      for (i = 0; i < m->vn[0]; i++)
        fprintf(f, "%.12g %.12g %.12g\n", i * m->c.length[0] / (m->vn[0] - 1),
                j * m->c.length[1] / (m->vn[1] - 1),
                k * m->c.length[2] / (m->vn[2] - 1));
  fprintf(f, "POINT_DATA %d\nVECTORS velocity_cm_year double\n", m->nv);
  for (k = 0; k < m->vn[2]; k++)
    for (j = 0; j < m->vn[1]; j++)
      for (i = 0; i < m->vn[0]; i++) {
        int id = node_id(i, j, k, m->vn);
        fprintf(f, "%.12g %.12g %.12g\n",
                m->velocity[3 * id] * m->c.velocity_scale,
                m->velocity[3 * id + 1] * m->c.velocity_scale,
                m->velocity[3 * id + 2] * m->c.velocity_scale);
      }
  fprintf(f, "SCALARS pressure_MPa double 1\nLOOKUP_TABLE default\n");
  for (k = 0; k < m->vn[2]; k++)
    for (j = 0; j < m->vn[1]; j++)
      for (i = 0; i < m->vn[0]; i++) {
        double x[3] = {i * m->h[0] / m->c.degree, j * m->h[1] / m->c.degree,
                       k * m->h[2] / m->c.degree},
               u[3], G[3][3], p;
        evaluate(m, x, u, G, &p);
        fprintf(f, "%.12g\n", p * m->P0 / 1e6);
      }
  fclose(f);
}
