
test_that("PL is correct for one persistence pair.", {
  # desired exact representation
  pl_exact <- list(cbind(c(-Inf, 0, 1, 2, Inf), c(0, 0, 1, 0, 0)))
  
  # double input
  pd <- matrix(c(0, 2), nrow = 1, ncol = 2)
  pl <- new(PersistenceLandscape, pd, TRUE, 0, 5, .001)
  expect_equal(pl$getInternal(), pl_exact)
  
  # integer input
  pd <- matrix(c(0L, 2L), nrow = 1, ncol = 2)
  pl <- new(PersistenceLandscape, pd, TRUE, 0, 5, .001)
  expect_equal(pl$getInternal(), pl_exact)
  
})
