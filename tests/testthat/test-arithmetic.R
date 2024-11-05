
# tiny example
pd <- cbind(start = 2, end = 5)
pl <- pl_new(pd, exact = TRUE)

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
})
