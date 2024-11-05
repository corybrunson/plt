
# tiny example
pd <- cbind(start = 0, end = 2)
pl <- pl_new(pd, exact = TRUE)

test_that("singletons are handled as lists", {
  # dist
  expect_no_error(pl_dist(pl))
})

test_that("norm agrees with integral", {
  # 1-norm
  expect_equal(pl_integrate(pl), pl_norm(pl, p = 1))
  # 2-norm
  expect_equal(pl_integrate(pl, p = 2) ^ (1/2), pl_norm(pl, p = 2))
  # 3-norm
  expect_equal(pl_integrate(pl, p = 3) ^ (1/3), pl_norm(pl, p = 3))
  # inf-norm?
})

# tinier example
pd2 <- cbind(start = 1, end = 2)
pl2 <- pl_new(pd2, exact = TRUE)

test_that("distance agrees with integral of difference", {
  # 1-norm
  expect_equal(pl_integrate(pl - pl2), pl_distance(pl, pl2, p = 1))
  # 2-norm
  expect_equal(pl_integrate(pl - pl2, p = 2)^(1/2), pl_distance(pl, pl2, p = 2))
  # 3-norm
  expect_equal(pl_integrate(pl - pl2, p = 3)^(1/3), pl_distance(pl, pl2, p = 3))
  # inf-norm?
})

# TODO: example with non-nested levels
