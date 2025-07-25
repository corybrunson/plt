set.seed(120246L)
t <- runif(12, 0, 2*pi)
x <- cbind(cos(t), sin(t))
pd_a <- TDA::alphaComplexDiag(x, maxdimension = 1)
pd_a$diagram[, c(2, 3)] <- sqrt(pd_a$diagram[, c(2, 3)])
pd_r <- ripserr::vietoris_rips(x, dim = 1)
pd_3 <- unclass(pd_a$diagram)
pd_2 <- as.matrix(pd_r)[, c(2, 3)]

test_that("`pl_new()` parses 2- and 3-column matrices correctly", {
  # 3-column, missing degree
  expect_error(pl_new(pd_3), regexp = "degree")
  # 3-column, degree specified
  expect_no_error(pl_new(pd_3, degree = 1))
  # 2-column, degree specified
  expect_warning(pl_new(pd_2, degree = 1), regexp = "degree")
  # 2-column, missing degree
  expect_no_error(pl_new(pd_2))
})

test_that("`pl_new()` accepts {TDA} output", {
  # full output
  expect_no_error(pl_new(pd_a, degree = 1))
  # diagram element
  expect_no_error(pl_new(pd_a$diagram, degree = 1))
})

test_that("`pl_new()` accepts {ripserr} output", {
  expect_no_error(pl_new(pd_r, degree = 1))
})

test_that("`pl_type()` works", {
  expect_no_error(pl_type(pl_new(pd_a,degree = 1)))
  expect_equal(pl_type(pl_new(pd_a,degree = 1)),
               pl_type(pl_new(pd_r,degree = 1)))
})

test_that("`pl_support()` works", {
  expect_no_error(pl_support(pl_new(pd_a,degree = 1)))
  expect_no_error(pl_support(pl_new(pd_r,degree = 1)))
})

test_that("`pl_is_exact()` works", {
  expect_equal(pl_is_exact(pl_new(pd_a,degree = 1)),
               pl_is_exact(pl_new(pd_r,degree = 1)))
})

test_that("`pl_limits()` works", {
  expect_no_error(pl_limits(pl_new(pd_a,degree = 1)))
  expect_no_error(pl_limits(pl_new(pd_r,degree = 1)))
  expect_equal(pl_limits(pl_new(pd_r,degree = 1))[2],
               pl_support(pl_new(pd_r,degree = 1))[2])
})

test_that("`pl_num_levels` can switch between exact and discrete PLs", {
  pd <- data.frame(dimension = 2, birth = c(1, 2, 3, 3), death = c(4, 4, 4, 6))
  pl <- pl_new(pd, degree = 2, exact = TRUE)
  pl2 <- pl_discretize(pl)
  
  expect_equal(pl_num_levels(pl), 4)
  expect_equal(pl_num_levels(pl2), 4)
  expect_error(pl_num_levels(pd))
  
})
