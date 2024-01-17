test_that("PL is correct for one persistence pair.", {
  # double input
  pd <- matrix(c(0, 2), nrow = 1, ncol = 2)
  pl <- new(PersistenceLandscape, pd, TRUE, 0, 5, .001)
  # integer input
  pd <- matrix(c(0L, 2L), nrow = 1, ncol = 2)
  pl <- new(PersistenceLandscape, pd, TRUE, 0, 5, .001)
})
