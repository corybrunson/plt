
pd <- matrix(c(0, 1, 0.5, 2, 0.75, 3), ncol = 2, byrow = TRUE)
pd2 <- matrix(c(0, 2), ncol = 2, byrow = TRUE)
pl_s <- pl_new(pd2, exact = FALSE)
pl <- pl_new(pd, exact = FALSE)


test_that("`pl_t` extracts the grid correctly", {
  t_vals <- pl_t(pl)
  expect_type(t_vals, "double")
  expect_true(is.numeric(t_vals))
  expect_true(length(t_vals) > 0)
  expect_equal(c(t_vals[1], t_vals[length(t_vals)]), pl_support(pl))
})

test_that("`pl_to_vector` returns the correct format and length", {
  vec <- pl_to_vector(pl, num_levels = 3)
  expect_type(vec, "double")
  expect_true(!is.null(attr(vec, "t_vals")))
  expect_equal(length(vec), 3 * length(attr(vec, "t_vals")))
})

test_that("`pl_to_vector` pads or trims levels correctly", {
  vec_short <- pl_to_vector(pl, num_levels = 2)
  vec_long <- pl_to_vector(pl, num_levels = 5)
  expect_equal(length(vec_short), 2 * length(attr(vec_short, "t_vals")))
  expect_equal(length(vec_long), 5 * length(attr(vec_long, "t_vals")))
  expect_equal(sum(vec_long[903:1505]), 0)
})

test_that("`pl_to_matrix` returns a matrix with correct dims", {
  mat <- pl_to_matrix(pl, pl_s)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 2)
  expect_equal(ncol(mat), max(length(pl_to_vector(pl_s)), length(pl_to_vector(pl))))
  expect_true(!is.null(attr(mat, "t_vals")))
})

test_that("`pl_from_vector` reconstructs a PL from a vector", {
  vec <- pl_to_vector(pl)
  pl2 <- pl_from_vector(vec)
  expect_true(inherits(pl2, "Rcpp_PersistenceLandscape"))
  expect_equal(attr(vec, "t_vals"), pl_t(pl2))
})

test_that("`pl_from_vector` respects drop_levels", {
  vec <- pl_to_vector(pl, num_levels = 5)
  pl2 <- pl_from_vector(vec, drop_levels = TRUE)
  expect_true(inherits(pl2, "Rcpp_PersistenceLandscape"))
  expect_lt(pl_num_levels(pl2), 5)
})

test_that("`pl_from_matrix` reconstructs multiple PLs correctly", {
  mat <- pl_to_matrix(pl, pl_s)
  pls <- pl_from_matrix(mat)
  expect_true(is.list(pls))
  expect_equal(length(pls), 2)
  expect_true(all(vapply(pls, inherits, TRUE, "Rcpp_PersistenceLandscape")))
  expect_equal(pl_distance(pls[[1]], pl), 0)
  expect_equal(pl_distance(pls[[2]], pl_s), 0)
})

test_that("`pl_from_matrix` handles single-row matrix as single PL", {
  mat <- pl_to_matrix(pl)
  pl_recon <- pl_from_matrix(mat)
  expect_true(inherits(pl_recon, "Rcpp_PersistenceLandscape"))
})