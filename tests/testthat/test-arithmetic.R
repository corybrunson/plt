
# tiny example
pd <- cbind(start = 2, end = 5)
pl <- pl_new(pd, exact = TRUE)
#larger example
pd2 <- matrix(c(0, 1, .25, 1.25), ncol = 2, byrow = TRUE)
pl2 <- pl_new(pd2, exact = FALSE)

test_that("singletons are handled as lists", {
  # sum
  expect_no_error(pl_sum(pl))
  # diff
  expect_no_error(pl_diff(pl))
  # mean
  expect_no_error(pl_mean(pl))
  # var
  expect_no_error(pl_var(pl))
  # sd
  expect_no_error(pl_sd(pl))
  # scale
  expect_error(pl_scale(pl))
  expect_no_error(pl_scale(pl,mult = 0))
  
})

test_that("classic functions work on PLs",{
  # max
  expect_no_error(pl_max(pl))
  # min
  expect_no_error(pl_min(pl))
  # moment
  expect_no_error(pl_moment(pl))
  # range
  expect_no_error(pl_range(pl))
})



test_that("vectorizations work",{
  # vmax
  expect_no_error(pl_vmax(pl))
  # vmin
  expect_no_error(pl_vmin(pl))
  # vmoment
  expect_no_error(pl_vmoment(pl))
  # vrange
  expect_no_error(pl_vrange(pl))
})

test_that("`pl_abs` correctly returns the absolute value of a PL", {
  pl3 <- -1 * pl
  pl4 <- pl_abs(pl3)
  
  expect_equal(pl_distance(pl, pl4), 0)
})

test_that("`pl_add` behaves as expected", {
  pla <- pl2 + pl2
  plb <- pl_add(pl2, pl2)
  plc <- 2 * pl2

  pls <- list(pla, plb, plc)
  expect_equal(sum(pl_dist(pls)), 0)
})

test_that("`pl_inner` obeys: symmetrical, nonnegative, and linearity ", {
  pl1 <- -1 * pl
  pla <- .75 * pl2
  
  expect_equal(pl_inner(pl1, pl2), pl_inner(pl2, pl1))
  expect_gte(pl_inner(pl2, pl2), 0)
  expect_equal(.75 * pl_inner(pl2, pl1), pl_inner(pla, pl1))
})
