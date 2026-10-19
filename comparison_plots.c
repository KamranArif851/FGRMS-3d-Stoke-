/* Only analytical-versus-iterative plots; no LU data are read. */
#include "plots.c"
static void convergence_plot(const char *prefix, int nodes) {
  char filename[2048], name[100], title[160], linebuf[512];
  snprintf(filename, sizeof(filename), "%s_history.csv", prefix);
  FILE *in = fopen(filename, "r");
  if (!in)
    fail("Missing convergence history");
  fgets(linebuf, sizeof(linebuf), in);
  double *est = allocate(200002, sizeof(double)),
         *truth = allocate(200002, sizeof(double));
  int *iter = allocate(200002, sizeof(int)), count = 0, maxiter = 1;
  double low = 0;
  while (count < 200002 && fgets(linebuf, sizeof(linebuf), in)) {
    char *end, *cursor = linebuf;
    iter[count] = (int)strtol(cursor, &end, 10);
    if (*end != ',')
      fail("Invalid history iteration");
    cursor = end + 1;
    est[count] = strtod(cursor, &end);
    if (*end != ',')
      fail("Invalid history estimate");
    cursor = end + 1;
    truth[count] = strncmp(cursor, "nan", 3) == 0 ? NAN : strtod(cursor, &end);
    low = fmin(low, log10(fmax(est[count], 1e-30)));
    maxiter = iter[count];
    count++;
  }
  if (!feof(in))
    fail("History exceeds supported iteration count");
  fclose(in);
  low = floor(low) - .3;
  if (maxiter < 1)
    maxiter = 1;
  snprintf(name, sizeof(name), "convergence_%d_nodes", nodes);
  snprintf(title, sizeof(title), "FGMRES convergence: %d velocity nodes",
           nodes);
  FILE *f = figure("figures", name, title, 1150, 650);
  axes(f, 120, 135, 940, 400, 0, maxiter, low, 0, "FGMRES iteration",
       "log10(relative residual)");
  for (int i = 1; i < count; i++)
    line(f, 120 + 940.0 * iter[i - 1] / maxiter,
         135 + 400 * log10(fmax(est[i - 1], 1e-30)) / low,
         120 + 940.0 * iter[i] / maxiter,
         135 + 400 * log10(fmax(est[i], 1e-30)) / low, red, 1.7);
  for (int i = 0; i < count; i++)
    if (isfinite(truth[i]))
      dot(f, 120 + 940.0 * iter[i] / maxiter,
          135 + 400 * log10(fmax(truth[i], 1e-30)) / low, 4, blue);
  text_svg(f, 36, 615, 16,
           "Red: Arnoldi estimate. Blue dots: recomputed true residual at "
           "restart boundaries and termination.",
           "black");
  finish(f, name, title);
  free(est);
  free(truth);
  free(iter);
}

int main(int argc, char **argv) {
  int is_study = argc == 1;
  if (argc != 1 && argc != 3)
    fail("Usage: comparison_plots.exe [input.cfg results/prefix]");
  const char *cfgs[] = {"results/q4_n2_input.cfg", "results/q4_n3_input.cfg",
                        "results/q4_443_input.cfg"};
  const char *pres[] = {"results/q4_n2_iter", "results/q4_n3_iter",
                        "results/q4_443_iter"};
  Model models[3];
  char prefix[2048];
  if (is_study)
    for (int i = 0; i < 3; i++)
      load_model(models + i, cfgs[i], pres[i]);
  else {
    snprintf(prefix, sizeof(prefix), "%s_iter", argv[2]);
    load_model(models, argv[1], prefix);
  }
  Model *m = models + (is_study ? 2 : 0);
  _mkdir("figures");
  gallery_folder = "figures";
  gallery = fopen("gallery.html", "w");
  if (!gallery)
    fail("Cannot open gallery");
  fprintf(gallery,
          "<!doctype html><meta charset='utf-8'><title>Analytical and FGMRES "
          "Stokes solution</title><style>body{font:17px "
          "Arial;max-width:1400px;margin:30px "
          "auto;padding:24px;color:#182942}p{line-height:1.6}.figures{display:"
          "grid;grid-template-columns:repeat(auto-fit,minmax(450px,1fr));gap:"
          "20px}figure{margin:0;padding:12px;border:1px solid "
          "#ddd}img{width:100%%}figcaption{padding:10px}</style><h1>Analytical "
          "solution and FGMRES</h1><p>Blue lines: analytical. Red dots: "
          "iterative numerical solution. Block Schur preconditioning and "
          "algebraic multigrid; no LU solve or comparison.</p><p>Domain: %.6g "
          "x %.6g x %.6g km. Elements: %d x %d x %d. Velocity degree: %d. "
          "Velocity nodes: %d; pressure nodes: %d.</p>",
          m->c.length[0], m->c.length[1], m->c.length[2], m->c.n[0], m->c.n[1],
          m->c.n[2], m->c.degree, m->nv, m->np);
  const char *last = is_study ? pres[2] : prefix;
  fprintf(
      gallery,
      "<p>Integrated relative L2 errors: velocity %.7g%%; pressure %.7g%%. "
      "These are FEM errors against the analytical solution, distinct from the "
      "algebraic solver residual.</p><p>For input instructions read "
      "START_HERE.txt. Every curve uses physical units as labelled. The 3-D "
      "surfaces are field heights over a middle-plane slice; streamlines show "
      "3-D flow.</p><h2>Convergence</h2><div class='figures'>",
      100 * stat(last, "velocity_relative_L2"),
      100 * stat(last, "pressure_relative_L2"));
  if (is_study)
    for (int i = 0; i < 3; i++)
      convergence_plot(pres[i], models[i].nv);
  else
    convergence_plot(prefix, m->nv);
  fprintf(gallery, "</div><h2>Velocity, pressure and stress fields</h2><div "
                   "class='figures'>");
  stress_summary(m, "figures");
  streamlines(m, "figures");
  surfaces(m, "figures");
  fprintf(gallery, "</div><h2>Analytical / iterative profiles and "
                   "errors</h2><div class='figures'>");
  profiles(m, "figures");
  if (is_study) {
    fprintf(gallery, "</div><h2>Mesh refinement study</h2><p>Fixed degree 4; "
                     "2x2x2, 3x3x3, 4x4x3 elements. The last refinement is "
                     "anisotropic.</p><div class='figures'>");
    study(models, pres, "figures");
  }
  fprintf(gallery, "</div>");
  fclose(gallery);
  for (int i = 0; i < (is_study ? 3 : 1); i++) {
    free(models[i].velocity);
    free(models[i].pressure);
  }
  puts("Saved analytical-versus-FGMRES figures and gallery.html");
  return 0;
}
