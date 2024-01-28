
# two sets of landscapes from similar but distinctive point clouds
set.seed(711018L)
circlescapes <- replicate(
  6,
  tdaunif::sample_circle(n = rpois(n = 1, lambda = 24)) |>
    ripserr::vietoris_rips(dim = 2L, threshold = 2) |>
    landscape(degree = 1, exact = TRUE)
)
toruscapes <- replicate(
  6,
  tdaunif::sample_torus_tube(n = rpois(n = 1, lambda = 24)) |>
    ripserr::vietoris_rips(dim = 2L, threshold = 2) |>
    landscape(degree = 1, exact = TRUE)
)

# null hypothesis tests
pl_z_test(circlescapes, toruscapes)
pl_perm_test(circlescapes, toruscapes)

\dontrun{
# benchmark one- and two-step computation of indicator-based linear form
bench::mark(
  pl$indicator_form(f, 0, p = 1),
  pl$indicator(f, 0)$integral(p = 1)
)
bench::mark(
  pl$indicator_form(f, r = 1, p = 1),
  pl$indicator(f, r = 1)$integral(p = 1)
)
}
