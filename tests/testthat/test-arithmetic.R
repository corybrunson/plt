
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
  # range
  expect_no_error(pl_range(pl))
  # scale ?
  
  
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

test_that
