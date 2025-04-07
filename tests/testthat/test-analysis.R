
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
  # inf-norm
  expect_equal(pl_vmax(pl_abs(pl)), pl_norm(pl, p = Inf))
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
  # inf-norm
  expect_equal(pl_vmax(pl_abs(pl - pl2)), pl_distance(pl, pl2, p = Inf))
})

# TODO: example with non-nested levels

# "Complex" example
pd3 <- data.frame(dim = 1, birth = c(0, 1), death = c(2, 2))
pl3 <- pl_new(pd3, degree = 1)

# Non-nested example
pd4 <- data.frame(dim = 1, birth = c(0, 2), death = c(4, 6))
pl4 <- pl_new(pd4, degree = 1)

test_that("`pl_indicator` applies indicator function correctly", {
  f <- list(c(0, 1)) # chosen so it splits PL in half
  
  pl_i <- pl_indicator(pl, f)
  
  expect_equal(pl_integrate(pl), 2 * pl_integrate(pl_i)) 
})

test_that("`pl_indicator` witholds `r` when applying `pl_integrate`", {
  f <- list(c(0, 2), c(1, 2))
  pl_i <- pl_indicator(pl3, f, 2)
  
  expect_equal(pl_integrate(pl_i), pl_integrate(pl) + pl_integrate(pl2) * 1/4)
})

test_that("`pl_indicator_form` with r = 0 matches `pl_indicator`", {
  f <- list(c(0, 1))
  pl_if <- pl_indicator_form(pl, f, r = 0, p = 1)
  pl_i <- pl_indicator(pl, f)
  
  f2 <- list(c(0, 2), c(1, 2))
  pl_if2 <- pl_indicator_form(pl3, f, r = 0, p = 1)
  pl_i2 <- pl_indicator(pl3, f, 2)
  
  expect_equal(pl_if, pl_integrate(pl_i)) 
  expect_equal(pl_if2, pl_integrate(pl_i2)) 
  
})

test_that("`pl_indicator_form` correctly applies `r`", {
  f <- list(c(0, 2), c(1, 2))
  pl_if <- pl_indicator_form(pl3, f, r = 2, p = 1)
  
  expect_equal(pl_if, pl_integrate(pl) + pl_integrate(pl2) * 1/4)
})

test_that("pl_indicator_form scales correctly with p", {
  f <- list(c(1, 2))
  
  # should all equal 1
  result_p1 <- pl_indicator_form(pl2, supports, r = 0, p = 1)
  result_p2 <- pl_indicator_form(pl2, supports, r = 0, p = 2) # not sure how this is working?
  result_pinf <- pl_indicator_form(pl2, supports, r = 0, p = Inf)
  

  expect_equal(result_p1, result_p2, tolerance = 1e-6)
  expect_equal(result_p2, result_pinf, tolerance = 1e-6)
  

  f <- list(c(0, 2), c(1, 2))
  # Weird results at this section
  result_t_p1 <- pl_indicator_form(pl3, f, r = 0, p = 1)
  result_t_p2 <- pl_indicator_form(pl3, f, r = 0, p = 2) # funny?
  result_t_pinf <- pl_indicator_form(pl3, f, r = 0, p = Inf)
  
  # Check expected ordering of norms
  expect_equal(result_t_p1, 1.25)
  expect_equal(result_t_p2, 1.0625)
  eexpect_equal(result_t_p1, result_t_pinf)
})

#pl_dist

test_that("`pl_dist` returns 0 matrix when using the same persistence landscape", {
  pl_list <- list(pl, pl, pl)
  distm <- pl_dist(pl_list)
  
  expect_equal(sum(distm), 0)
}
)
