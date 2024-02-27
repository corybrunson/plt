
test_that("PL is correct for one persistence pair", {
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

test_that("delimiting a discrete PL yields expected grids", {
  pd <- matrix(c(0L, 2L), nrow = 1, ncol = 2)
  pl <- new(PersistenceLandscape, pd, FALSE, 0, 2.05, .1)
  
  # delimit to same max (with error)
  expect_equal(
    dim(pl_delimit(pl, xmax = 2.05)$getInternal())[[2L]],
    dim(pl$getInternal())[[2L]]
  )
  expect_equal(
    dim(pl_delimit(pl, xmax = 2.04)$getInternal())[[2L]],
    dim(pl$getInternal())[[2L]]
  )
  expect_equal(
    dim(pl_delimit(pl, xmax = 2.06)$getInternal())[[2L]],
    dim(pl$getInternal())[[2L]]
  )
  
  # delimit to lower max
  expect_equal(
    dim(pl_delimit(pl, xmax = 2)$getInternal())[[2L]],
    dim(pl$getInternal())[[2L]] - 1L
  )
  # error: new max excludes some support
  expect_error(
    dim(pl_delimit(pl, xmax = 1.99)$getInternal())[[2L]],
    "support"
  )

  # delimit to higher max
  expect_equal(
    dim(pl_delimit(pl, xmax = 2.15)$getInternal())[[2L]],
    dim(pl$getInternal())[[2L]] + 1L
  )
  expect_equal(
    dim(pl_delimit(pl, xmax = 2.11)$getInternal())[[2L]],
    dim(pl$getInternal())[[2L]] + 1L
  )
})
