/* All figure geometry, exact fields, errors and SVG output are computed in C.
 * SVG is a standard vector image; open gallery.html in a browser to view it.
 * No Python, PyTorch, JavaScript, plotting package, or network is needed. */
#define main stokes_solver_entry
#include "stokes_fem.c"
#undef main
#include <direct.h>

static const char *fields[17] = {
    "vx",         "vy",         "vz",         "speed",      "pressure",
    "viscous_xx", "viscous_yy", "viscous_zz", "viscous_xy", "viscous_xz",
    "viscous_yz", "total_xx",   "total_yy",   "total_zz",   "total_xy",
    "total_xz",   "total_yz"};
static const char *units(int f) { return f < 4 ? "cm/year" : "MPa"; }
static const char *notation[17] = {
    "v_x",  "v_y",  "v_z", "|v|", "P",   "σ′xx", "σ′yy", "σ′zz", "σ′xy",
    "σ′xz", "σ′yz", "σxx", "σyy", "σzz", "σxy",  "σxz",  "σyz"};
static const char *blue = "#1459dc", *red = "#df253b";
static FILE *gallery;
static const char *gallery_folder = "figures";
static void text_svg(FILE *f, double x, double y, int sz, const char *s,
                     const char *color) {
  fprintf(f, "<text x='%.2f' y='%.2f' font-size='%d' fill='%s'>", x, y, sz,
          color);
  while (*s) {
    if (strncmp(s, "\xCF\x83", 2) == 0) {
      fputs("\xCF\x83", f);
      s += 2;
      if (strncmp(s, "\xE2\x80\xB2", 3) == 0) {
        fputs("\xE2\x80\xB2", f);
        s += 3;
      }
      if (s[0] && s[1] && strchr("xyz", s[0]) && strchr("xyz", s[1])) {
        fprintf(f, "<tspan baseline-shift='sub' font-size='70%%'>%c%c</tspan>",
                s[0], s[1]);
        s += 2;
      }
    } else
      fputc(*s++, f);
  }
  fputs("</text>\n", f);
}
static FILE *figure(const char *dir, const char *name, const char *title, int w,
                    int h) {
  char path[2048];
  snprintf(path, sizeof(path), "%s/%s.svg", dir, name);
  FILE *f = fopen(path, "w");
  if (!f)
    fail("Cannot open SVG");
  fprintf(f,
          "<svg xmlns='http://www.w3.org/2000/svg' width='%d' height='%d' "
          "viewBox='0 0 %d %d'><rect width='100%%' height='100%%' "
          "fill='white'/><g font-family='Times New Roman,serif'>\n",
          w, h, w, h);
  text_svg(f, 36, 38, 24, title, "#16233c");
  return f;
}
static void finish(FILE *f, const char *name, const char *title) {
  fprintf(f, "</g></svg>");
  fclose(f);
  if (gallery)
    fprintf(gallery,
            "<figure><a href='%s/%s.svg'><img loading='lazy' "
            "src='%s/%s.svg' "
            "alt='%s'></a><figcaption>%s</figcaption></figure>\n",
            gallery_folder, name, gallery_folder, name, title, title);
}
static double stat(const char *prefix, const char *key) {
  char path[2048], buf[16384], needle[128];
  snprintf(path, sizeof(path), "%s_stats.json", prefix);
  FILE *f = fopen(path, "r");
  if (!f)
    fail("Missing statistics file");
  size_t n = fread(buf, 1, sizeof(buf) - 1, f);
  buf[n] = 0;
  fclose(f);
  snprintf(needle, sizeof(needle), "\"%s\":", key);
  char *p = strstr(buf, needle);
  if (!p)
    fail("Missing JSON metric");
  return strtod(p + strlen(needle), NULL);
}
static void load_model(Model *m, const char *cfg, const char *prefix) {
  char status_path[2048],status_text[2048];
  size_t prefix_length=strlen(prefix);
  if(prefix_length<5 || strcmp(prefix+prefix_length-5,"_iter"))fail("Expected iterative results prefix ending _iter");
  snprintf(status_path,sizeof(status_path),"%.*s_solver.json",(int)(prefix_length-5),prefix);
  FILE *status=fopen(status_path,"r");if(!status)fail("Missing solver status");
  size_t status_size=fread(status_text,1,sizeof(status_text)-1,status);status_text[status_size]=0;fclose(status);
  if(!strstr(status_text,"\"converged\":true"))fail("The latest solve did not complete successfully; do not plot stale fields");
  memset(m, 0, sizeof(*m));
  m->c = defaults();
  read_config(&m->c, cfg);
  m->L0 = fmax(m->c.length[0], fmax(m->c.length[1], m->c.length[2])) * 1000;
  m->U0 = m->c.velocity_scale * .01 / YEAR;
  m->P0 = m->c.eta0 * m->U0 / m->L0;
  m->F0 = m->P0 / m->L0;
  m->nv = m->np = 1;
  for (int d = 0; d < 3; d++) {
    m->ell[d] = m->c.length[d] * 1000 / m->L0;
    m->h[d] = m->ell[d] / m->c.n[d];
    m->vn[d] = m->c.degree * m->c.n[d] + 1;
    m->pn[d] = (m->c.degree - 1) * m->c.n[d] + 1;
    m->nv *= m->vn[d];
    m->np *= m->pn[d];
  }
  m->local_v = (m->c.degree + 1) * (m->c.degree + 1) * (m->c.degree + 1);
  m->local_p = m->c.degree * m->c.degree * m->c.degree;
  m->local_dof = 3 * m->local_v;
  m->velocity = allocate(3 * m->nv, sizeof(double));
  m->pressure = allocate(m->np, sizeof(double));
  char path[2048], line[1024];
  double x, y, z;
  snprintf(path, sizeof(path), "%s_velocity.csv", prefix);
  FILE *f = fopen(path, "r");
  if (!f)
    fail("Missing velocity CSV");
  fgets(line, sizeof(line), f);
  for (int i = 0; i < m->nv; i++)
    if (fscanf(f, "%lf,%lf,%lf,%lf,%lf,%lf", &x, &y, &z, &m->velocity[3 * i],
               &m->velocity[3 * i + 1], &m->velocity[3 * i + 2]) != 6)
      fail("Invalid velocity CSV");
  fclose(f);
  snprintf(path, sizeof(path), "%s_pressure.csv", prefix);
  f = fopen(path, "r");
  if (!f)
    fail("Missing pressure CSV");
  fgets(line, sizeof(line), f);
  for (int i = 0; i < m->np; i++)
    if (fscanf(f, "%lf,%lf,%lf,%lf", &x, &y, &z, &m->pressure[i]) != 4)
      fail("Invalid pressure CSV");
  fclose(f);
  if (!has_exact(m))
    fail(
        "These comparison figures require a valid manufactured exact solution");
}
static void sample(const Model *m, double X, double Y, double Z, double *v,
                   double *e) {
  double x[3] = {X * m->ell[0], Y * m->ell[1], Z * m->ell[2]};
  all_fields(m, x, v, e);
}
static void line(FILE *f, double x1, double y1, double x2, double y2,
                 const char *color, double width) {
  fprintf(f,
          "<path d='M%.2f %.2fL%.2f %.2f' fill='none' stroke='%s' "
          "stroke-width='%.2f'/>\n",
          x1, y1, x2, y2, color, width);
}
static void dot(FILE *f, double x, double y, double r, const char *color) {
  fprintf(f, "<circle cx='%.2f' cy='%.2f' r='%.2f' fill='%s'/>\n", x, y, r,
          color);
}
static void axes(FILE *f, double l, double t, double w, double h, double xmin,
                 double xmax, double ymin, double ymax, const char *xl,
                 const char *yl) {
  line(f, l, t + h, l + w, t + h, "#607087", 1);
  line(f, l, t, l, t + h, "#607087", 1);
  char s[80];
  for (int i = 0; i < 5; i++) {
    double r = i / 4.0;
    snprintf(s, sizeof(s), "%.3g", xmin + r * (xmax - xmin));
    text_svg(f, l + r * w - 12, t + h + 23, 12, s, "#35465f");
    snprintf(s, sizeof(s), "%.3g", ymin + r * (ymax - ymin));
    text_svg(f, l - 65, t + h - r * h + 4, 12, s, "#35465f");
  }
  text_svg(f, l + w / 2 - 20, t + h + 47, 14, xl, "#35465f");
  text_svg(f, l, t - 13, 14, yl, "#35465f");
}
static void curve(FILE *f, const double *y, int n, double l, double t, double w,
                  double h, double lo, double hi, const char *color, int dots) {
  if (!dots)
    fprintf(f, "<path d='");
  for (int i = 0; i < n; i++) {
    double x = l + w * i / (n - 1.0), v = t + h * (hi - y[i]) / (hi - lo);
    if (dots) {
      if (i % 5 == 0)
        dot(f, x, v, 2.6, color);
    } else
      fprintf(f, "%c%.2f %.2f", i ? 'L' : 'M', x, v);
  }
  if (!dots)
    fprintf(f,
            "' fill='none' stroke='%s' stroke-width='2.3' "
            "stroke-linejoin='round'/>\n",
            color);
}
static void profiles(Model *m, const char *dir) {
  double v[3][201][17], e[3][201][17];
  for (int a = 0; a < 3; a++)
    for (int i = 0; i <= 200; i++) {
      double x[3] = {.31, .37, .43};
      x[a] = i / 200.0;
      sample(m, x[0], x[1], x[2], v[a][i], e[a][i]);
    }
  for (int k = 0; k < 17; k++)
    for (int a = 0; a < 3; a++) {
      char name[128], title[160], unit[80];
      snprintf(name, sizeof(name), "%s_%c_profile", fields[k], 'x' + a);
      snprintf(title, sizeof(title), "%s | %c-direction profile and error",
               notation[k], 'x' + a);
      FILE *f = figure(dir, name, title, 1100, 780);
      text_svg(
          f, 36, 69, 14,
          "Exact: blue line     FGMRES: red dots     No mesh/grid lines",
          "#35465f");
      double n[201], ex[201], err[201], lo = 1e300, hi = -1e300, emax = 0;
      for (int i = 0; i <= 200; i++) {
        n[i] = v[a][i][k];
        ex[i] = e[a][i][k];
        err[i] = n[i] - ex[i];
        lo = fmin(lo, fmin(n[i], ex[i]));
        hi = fmax(hi, fmax(n[i], ex[i]));
        emax = fmax(emax, fabs(err[i]));
      }
      double pad = fmax((hi - lo) * .08, 1e-12);
      lo -= pad;
      hi += pad;
      emax = fmax(emax * 1.08, 1e-16);
      snprintf(unit, sizeof(unit), "%s (%s)", notation[k], units(k));
      char xl[40];
      snprintf(xl, sizeof(xl), "%c (km)", 'x' + a);
      axes(f, 110, 120, 920, 320, 0, m->c.length[a], lo, hi, xl, unit);
      curve(f, ex, 201, 110, 120, 920, 320, lo, hi, blue, 0);
      curve(f, n, 201, 110, 120, 920, 320, lo, hi, red, 1);
      snprintf(unit, sizeof(unit), "Numerical minus exact (%s)", units(k));
      axes(f, 110, 550, 920, 145, 0, m->c.length[a], -emax, emax, xl, unit);
      line(f, 110, 622.5, 1030, 622.5, "#ccd4e0", 1);
      curve(f, err, 201, 110, 550, 920, 145, -emax, emax, red, 0);
      text_svg(f, 36, 765, 12,
               "Fixed coordinates (where applicable): x/Lx = 0.31, y/Ly = "
               "0.37, z/Lz = 0.43.",
               "#607087");
      finish(f, name, title);
    }
}
typedef struct {
  double x, y;
} Point;
static Point project(double x, double y, double z, double cx, double cy,
                     double scale) {
  Point p = {cx + scale * .866 * (x - y),
             cy + scale * (.37 * (x + y) - .95 * z)};
  return p;
}
static void color(double value, double shade, int palette, char out[16]) {
  double r, g, b, t = fmax(0, fmin(1, value));
  if (palette == 0) {
    r = 35 + 115 * t;
    g = 82 + 125 * t;
    b = 177 + 70 * t;
  } else if (palette == 1) {
    r = 164 + 75 * t;
    g = 40 + 116 * t;
    b = 66 + 104 * t;
  } else {
    r = 45 + 205 * t;
    g = 20 + 195 * t * t;
    b = 85 - 40 * t;
  }
  snprintf(out, 16, "#%02x%02x%02x", (int)fmin(255, r * shade),
           (int)fmin(255, g * shade), (int)fmin(255, b * shade));
}
static void surface(FILE *f, double *z, int N, double cx, double cy,
                    double scale, double lo, double hi, int palette,
                    const Model *m, const char *unit) {
  double span = fmax(hi - lo, 1e-30);
  for (int sum = 0; sum <= 2 * (N - 2); sum++)
    for (int i = 0; i < N - 1; i++) {
      int j = sum - i;
      if (j < 0 || j >= N - 1)
        continue;
      int ids[4] = {i * N + j, (i + 1) * N + j, (i + 1) * N + j + 1,
                    i * N + j + 1};
      double mean = 0;
      Point p[4];
      for (int d = 0; d < 4; d++) {
        int a = ids[d] / N, b = ids[d] % N;
        double value = (z[ids[d]] - lo) / span;
        mean += value / 4;
        p[d] =
            project(a / (N - 1.0), b / (N - 1.0), .95 * value, cx, cy, scale);
      }
      double dx = (z[ids[1]] - z[ids[0]]) / span * (N - 1),
             dy = (z[ids[3]] - z[ids[0]]) / span * (N - 1);
      double light = .72 + .25 * fmax(0, (-.35 * dx - .45 * dy + 1) /
                                             sqrt(1 + dx * dx + dy * dy));
      char col[16];
      color(mean, light, palette, col);
      fprintf(f,
              "<polygon points='%.2f,%.2f %.2f,%.2f %.2f,%.2f %.2f,%.2f' "
              "fill='%s' stroke='%s' stroke-width='.35'/>\n",
              p[0].x, p[0].y, p[1].x, p[1].y, p[2].x, p[2].y, p[3].x, p[3].y,
              col, col);
    }
  Point o = project(1, 1, 0, cx, cy, scale),
        a = project(0, 1, 0, cx, cy, scale),
        b = project(1, 0, 0, cx, cy, scale),
        top = project(0, 1, 1.05, cx, cy, scale);
  line(f, o.x, o.y, a.x, a.y, "black", 1);
  line(f, o.x, o.y, b.x, b.y, "black", 1);
  line(f, a.x, a.y, top.x, top.y, "black", 1);
  for (int tick = 0; tick < 3; tick++) {
    double r = tick / 2.0;
    char value[64];
    Point p = project(r, 1, 0, cx, cy, scale);
    line(f, p.x, p.y, p.x - 3, p.y + 5, "black", 1);
    snprintf(value, sizeof(value), "%.0f", r * m->c.length[0]);
    text_svg(f, p.x - 35, p.y + 20, 13, value, "black");
    p = project(1, r, 0, cx, cy, scale);
    line(f, p.x, p.y, p.x + 3, p.y + 5, "black", 1);
    snprintf(value, sizeof(value), "%.0f", r * m->c.length[1]);
    text_svg(f, p.x + 14, p.y + 17, 13, value, "black");
    p = project(0, 1, .95 * r, cx, cy, scale);
    line(f, p.x - 5, p.y, p.x, p.y, "black", 1);
    snprintf(value, sizeof(value), "%.3g", lo + (hi - lo) * r);
    text_svg(f, p.x - 63, p.y + 4, 13, value, "black");
  }
  text_svg(f, cx - 95, o.y + 28, 15, "x (km)", "black");
  text_svg(f, cx + 65, o.y + 28, 15, "y (km)", "black");
  text_svg(f, top.x - 35, top.y - 12, 15, unit, "black");
  char s[120];
  snprintf(s, sizeof(s), "Height range: %.4g to %.4g", lo, hi);
  text_svg(f, cx - 170, cy + scale * .9 + 30, 13, s, "#35465f");
}
static void surfaces(Model *m, const char *dir) {
  enum { N = 81 };
  double *nv = allocate(N * N * 17, sizeof(double)),
         *ev = allocate(N * N * 17, sizeof(double));
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      sample(m, i / (N - 1.0), j / (N - 1.0), .43, nv + (i * N + j) * 17,
             ev + (i * N + j) * 17);
  double *n = allocate(N * N, sizeof(double)),
         *e = allocate(N * N, sizeof(double)),
         *err = allocate(N * N, sizeof(double));
  for (int k = 0; k < 17; k++) {
    double lo = 1e300, hi = -1e300, emax = 0;
    for (int i = 0; i < N * N; i++) {
      n[i] = nv[i * 17 + k];
      e[i] = ev[i * 17 + k];
      err[i] = fabs(n[i] - e[i]);
      lo = fmin(lo, fmin(n[i], e[i]));
      hi = fmax(hi, fmax(n[i], e[i]));
      emax = fmax(emax, err[i]);
    }
    char name[128], title[160], caption[160];
    snprintf(name, sizeof(name), "%s_smooth_3d", fields[k]);
    snprintf(title, sizeof(title),
             "%s (%s) | exact, numerical and absolute error", notation[k],
             units(k));
    FILE *f = figure(dir, name, title, 1500, 620);
    text_svg(f, 36, 68, 14,
             "Slice at z/Lz = 0.43; vertical height represents field value, "
             "not physical z. Exact and numerical share a height scale.",
             "#607087");
    text_svg(f, 150, 110, 20, "(a) Exact solution", blue);
    text_svg(f, 650, 110, 20, "(b) FGMRES solution", red);
    text_svg(f, 1150, 110, 20, "(c) Absolute error", "#653b85");
    surface(f, e, N, 270, 340, 180, lo, hi, 0, m, units(k));
    surface(f, n, N, 770, 340, 180, lo, hi, 1, m, units(k));
    surface(f, err, N, 1270, 340, 180, 0, fmax(emax, 1e-30), 2, m, units(k));
    snprintf(
        caption, sizeof(caption),
        "Error height uses its own labelled scale. No smoothing of numerical "
        "stresses; element-face jumps remain visible. Units: %s.",
        units(k));
    text_svg(f, 36, 602, 13, caption, "#607087");
    finish(f, name, title);
  }
  free(n);
  free(e);
  free(err);
  free(nv);
  free(ev);
}
static int velocity_at(Model *m, const double r[3], double v[3],
                       int numerical) {
  for (int i = 0; i < 3; i++)
    if (r[i] < 0 || r[i] > 1)
      return 0;
  double x[3], G[3][3], p, lap[3], gp[3];
  for (int i = 0; i < 3; i++)
    x[i] = r[i] * m->ell[i];
  if (numerical)
    evaluate(m, x, v, G, &p);
  else
    exact(m, x, v, G, lap, &p, gp);
  for (int i = 0; i < 3; i++)
    v[i] /= m->ell[i];
  return 1;
}
static void streamlines(Model *m, const char *dir) {
  const char *name = "velocity_streamlines_3d",
             *title = "Velocity streamlines in the full 3-D domain";
  FILE *f = figure(dir, name, title, 1400, 800);
  text_svg(f, 36, 70, 14,
           "Same seeds and RK4 integration in both panels. Curves follow the "
           "velocity field; no mesh points or background grid.",
           "#607087");
  for (int numerical = 0; numerical < 2; numerical++) {
    double cx = 350 + 700 * numerical, cy = 430, scale = 310;
    text_svg(f, cx - 100, 118, 22,
             numerical ? "FGMRES solution" : "Exact solution",
             numerical ? red : blue);
    for (int sy = 0; sy < 4; sy++)
      for (int sz = 0; sz < 4; sz++) {
        double r[3] = {.42, .18 + .20 * sy, .18 + .20 * sz};
        fprintf(f, "<path d='");
        int count = 0;
        for (int step = 0; step < 420; step++) {
          double a[3], b[3], c[3], d[3], tmp[3], dt = .016;
          if (!velocity_at(m, r, a, numerical))
            break;
          Point p = project(r[0], r[1], r[2], cx, cy, scale);
          fprintf(f, "%c%.2f %.2f", count++ ? 'L' : 'M', p.x, p.y);
          for (int k = 0; k < 3; k++)
            tmp[k] = r[k] + dt * a[k] / 2;
          if (!velocity_at(m, tmp, b, numerical))
            break;
          for (int k = 0; k < 3; k++)
            tmp[k] = r[k] + dt * b[k] / 2;
          if (!velocity_at(m, tmp, c, numerical))
            break;
          for (int k = 0; k < 3; k++)
            tmp[k] = r[k] + dt * c[k];
          if (!velocity_at(m, tmp, d, numerical))
            break;
          for (int k = 0; k < 3; k++)
            r[k] += dt * (a[k] + 2 * b[k] + 2 * c[k] + d[k]) / 6;
        }
        fprintf(f,
                "' fill='none' stroke='%s' stroke-opacity='.72' "
                "stroke-width='2.2' stroke-linecap='round'/>\n",
                numerical ? red : blue);
      }
    Point o = project(0, 0, 0, cx, cy, scale);
    for (int k = 0; k < 3; k++) {
      double r[3] = {0};
      r[k] = 1.05;
      Point p = project(r[0], r[1], r[2], cx, cy, scale);
      line(f, o.x, o.y, p.x, p.y, "#6c809d", 1);
      char label[80];
      snprintf(label, sizeof(label), "%c: %.0f km", 'x' + k, m->c.length[k]);
      text_svg(f, p.x - 25, p.y + (k == 2 ? -10 : 25), 14, label, "#35465f");
    }
  }
  text_svg(f, 36, 770, 13,
           "Streamlines show direction and geometry. Quantitative velocity, "
           "speed and error values are in the profile and surface figures.",
           "#607087");
  finish(f, name, title);
}
static void study(Model *models, const char **prefix, const char *dir) {
  const char *metrics[] = {"velocity_relative_L2",
                           "pressure_relative_L2",
                           "velocity_L2",
                           "pressure_L2",
                           "velocity_H1_semi",
                           "speed_L2",
                           "strain_L2",
                           "viscous_stress_L2",
                           "total_stress_L2",
                           "energy_error",
                           "vx_L2",
                           "vy_L2",
                           "vz_L2",
                           "velocity_Linf_sampled",
                           "pressure_Linf_sampled",
                           "speed_Linf_sampled",
                           "divergence_L2",
                           "scaled_residual"};
  FILE *csv = fopen("node_comparison.csv", "w");
  if (!csv)
    fail("Cannot open study CSV");
  fprintf(csv, "velocity_nodes,pressure_nodes,nx,ny,nz");
  for (int k = 0; k < 18; k++)
    fprintf(csv, ",%s", metrics[k]);
  fprintf(csv, "\n");
  for (int j = 0; j < 3; j++) {
    fprintf(csv, "%d,%d,%d,%d,%d", models[j].nv, models[j].np, models[j].c.n[0],
            models[j].c.n[1], models[j].c.n[2]);
    for (int k = 0; k < 18; k++)
      fprintf(csv, ",%.17g", stat(prefix[j], metrics[k]));
    fprintf(csv, "\n");
  }
  fclose(csv);
  for (int k = 0; k < 18; k++) {
    double values[3], lo = 1e300, hi = -1e300;
    for (int j = 0; j < 3; j++) {
      values[j] =
          log10(fmax(stat(prefix[j], metrics[k]) * (k < 2 ? 100 : 1), 1e-30));
      lo = fmin(lo, values[j]);
      hi = fmax(hi, values[j]);
    }
    lo -= .15;
    hi += .15;
    char name[128], title[160];
    snprintf(name, sizeof(name), "nodes_%s", metrics[k]);
    snprintf(title, sizeof(title), "Node study | %s", metrics[k]);
    FILE *f = figure(dir, name, title, 1100, 650);
    text_svg(f, 36, 70, 14,
             "Fixed Q4/Q3 order. Mesh dimensions are stated in the report; "
             "largest mesh refines two directions.",
             "#607087");
    axes(f, 130, 140, 860, 390, models[0].nv, models[2].nv, lo, hi,
         "Velocity nodes",
         k < 2 ? "log10(relative error in %)"
               : "log10(nondimensional error / residual)");
    for (int j = 0; j < 3; j++) {
      double x = 130 + 860.0 * (models[j].nv - models[0].nv) /
                           (models[2].nv - models[0].nv),
             y = 140 + 390 * (hi - values[j]) / (hi - lo);
      if (j) {
        double xp = 130 + 860.0 * (models[j - 1].nv - models[0].nv) /
                              (models[2].nv - models[0].nv),
               yp = 140 + 390 * (hi - values[j - 1]) / (hi - lo);
        line(f, xp, yp, x, y, red, 2);
      }
      dot(f, x, y, 5, red);
      char s[80];
      snprintf(s, sizeof(s), "%.4g%s (%d nodes)", pow(10, values[j]),
               k < 2 ? "%" : "", models[j].nv);
      text_svg(f, fmin(x, 800), y - 15, 14, s, "#35465f");
    }
    finish(f, name, title);
  }
  for (int k = 0; k < 17; k++) {
    double n[3][201], e[201], lo = 1e300, hi = -1e300;
    for (int j = 0; j < 3; j++)
      for (int i = 0; i <= 200; i++) {
        double v[17], ex[17];
        sample(&models[j], i / 200.0, .37, .43, v, ex);
        n[j][i] = v[k];
        e[i] = ex[k];
        lo = fmin(lo, fmin(v[k], ex[k]));
        hi = fmax(hi, fmax(v[k], ex[k]));
      }
    double pad = fmax(.08 * (hi - lo), 1e-12);
    lo -= pad;
    hi += pad;
    char name[128], title[160];
    snprintf(name, sizeof(name), "nodes_%s_profiles", fields[k]);
    snprintf(title, sizeof(title), "%s | comparison at different node counts",
             notation[k]);
    FILE *f = figure(dir, name, title, 1500, 560);
    for (int j = 0; j < 3; j++) {
      char s[100];
      snprintf(s, sizeof(s), "%d velocity / %d pressure nodes", models[j].nv,
               models[j].np);
      text_svg(f, 100 + 500 * j, 85, 16, s, "#35465f");
      axes(f, 90 + 500 * j, 140, 380, 300, 0, models[j].c.length[0], lo, hi,
           "x (km)", units(k));
      curve(f, e, 201, 90 + 500 * j, 140, 380, 300, lo, hi, blue, 0);
      curve(f, n[j], 201, 90 + 500 * j, 140, 380, 300, lo, hi, red, 1);
    }
    text_svg(f, 36, 535, 14,
             "Exact: blue line. FGMRES: red dots. Identical vertical scales "
             "across panels.",
             "#35465f");
    finish(f, name, title);
  }
}
static void stress_summary(Model *m, const char *dir) {
  double numerical[201][17], analytical[201][17];
  for (int i = 0; i <= 200; i++)
    sample(m, i / 200.0, .37, .43, numerical[i], analytical[i]);
  for (int total = 0; total < 2; total++) {
    const char *name = total ? "six_total_stresses_and_pressure"
                             : "six_viscous_stresses_and_pressure";
    const char *title = total ? "Total stress tensor components and pressure"
                              : "Viscous stress tensor components and pressure";
    FILE *f = figure(dir, name, title, 1500, 1180);
    text_svg(f, 36, 70, 17,
             "Exact: blue line. FGMRES: red dots. All stresses and pressure "
             "are in MPa.",
             "black");
    for (int panel = 0; panel < 7; panel++) {
      int k = panel < 6 ? (total ? 11 : 5) + panel : 4;
      double n[201], e[201], lo = 1e300, hi = -1e300;
      for (int i = 0; i <= 200; i++) {
        n[i] = numerical[i][k];
        e[i] = analytical[i][k];
        lo = fmin(lo, fmin(n[i], e[i]));
        hi = fmax(hi, fmax(n[i], e[i]));
      }
      double padding = fmax(.1 * (hi - lo), 1e-12);
      lo -= padding;
      hi += padding;
      double left = 90 + 500 * (panel % 3), top = 155 + 340 * (panel / 3);
      char label[100];
      snprintf(label, sizeof(label), "(%c) %s", 'a' + panel, notation[k]);
      text_svg(f, left, top - 42, 22, label, "black");
      axes(f, left, top, 370, 220, 0, m->c.length[0], lo, hi, "x (km)", "MPa");
      curve(f, e, 201, left, top, 370, 220, lo, hi, blue, 0);
      curve(f, n, 201, left, top, 370, 220, lo, hi, red, 1);
    }
    text_svg(f, 610, 865, 24, "σ = σ′ − P I", "black");
    text_svg(f, 610, 915, 19,
             "σxx = σ′xx − P;   σyy = σ′yy − P;   σzz = σ′zz − P", "black");
    text_svg(f, 610, 955, 19, "σxy = σ′xy;   σxz = σ′xz;   σyz = σ′yz",
             "black");
    text_svg(f, 610, 1005, 17,
             "Symmetry gives six independent stress components.", "black");
    text_svg(f, 36, 1150, 16,
             "Profile at y/Ly = 0.37 and z/Lz = 0.43. Component signs and "
             "pressure gauge are unchanged.",
             "black");
    finish(f, name, title);
  }
  FILE *f =
      figure(dir, "stress_tensor_and_pressure",
             "Stress notation, pressure and the Stokes equation", 1400, 1000);
  text_svg(f, 55, 95, 23, "div(σ′) − grad(P) + f = 0;    div(v) = 0", "black");
  text_svg(f, 55, 145, 23, "σ′ij = η (∂vi/∂xj + ∂vj/∂xi);       σ = σ′ − P I",
           "black");
  int tensor[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
  for (int total = 0; total < 2; total++) {
    double left = 100 + 650 * total;
    text_svg(f, left, 205, 24, total ? "Total stress σ" : "Viscous stress σ′",
             "black");
    line(f, left - 20, 230, left - 20, 400, "black", 2);
    line(f, left + 440, 230, left + 440, 400, "black", 2);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        text_svg(f, left + 145 * j, 270 + 55 * i, 25,
                 notation[(total ? 11 : 5) + tensor[i][j]], "black");
  }
  text_svg(f, 55, 450, 18,
           "Pressure subtracts from the three diagonal entries only. Shear "
           "entries are unchanged.",
           "black");
  text_svg(f, 55, 495, 21,
           "Values at x/Lx = 0.31, y/Ly = 0.37, z/Lz = 0.43 (MPa)", "black");
  text_svg(f, 60, 540, 18, "Quantity", "black");
  text_svg(f, 280, 540, 18, "Exact", blue);
  text_svg(f, 520, 540, 18, "Numerical", red);
  text_svg(f, 800, 540, 18, "Numerical − exact", "black");
  double v[17], e[17];
  sample(m, .31, .37, .43, v, e);
  FILE *csv = fopen("stress_probe.csv", "w");
  if (!csv)
    fail("Cannot write stress probe");
  fprintf(csv, "quantity,exact_MPa,numerical_MPa,signed_error_MPa\n");
  for (int row = 0; row < 7; row++) {
    int k = row < 6 ? 5 + row : 4;
    char a[80], b[80], c[80];
    snprintf(a, sizeof(a), "%.9g", e[k]);
    snprintf(b, sizeof(b), "%.9g", v[k]);
    snprintf(c, sizeof(c), "%.5g", v[k] - e[k]);
    text_svg(f, 60, 585 + 48 * row, 20, notation[k], "black");
    text_svg(f, 280, 585 + 48 * row, 20, a, blue);
    text_svg(f, 520, 585 + 48 * row, 20, b, red);
    text_svg(f, 800, 585 + 48 * row, 20, c, "black");
    fprintf(csv, "%s,%.17g,%.17g,%.17g\n", fields[k], e[k], v[k], v[k] - e[k]);
  }
  fclose(csv);
  text_svg(f, 55, 970, 16,
           "Separate total-stress profiles, all three profile directions, 3-D "
           "surfaces and signed-error panels are also supplied.",
           "black");
  finish(f, "stress_tensor_and_pressure",
         "Stress notation, pressure and the Stokes equation");
}
