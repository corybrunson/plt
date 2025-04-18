
test_that("Addition expands domain as needed.", {
  
  # exact landscapes
  pl1 <- pl_new(matrix(c(1, 3), ncol = 2L), exact = TRUE)
  pl2 <- pl_new(matrix(c(2, 4), ncol = 2L), exact = TRUE, xby = 0.2)
  pl12 <- pl1 + pl2
  expect_equal(pl12$xMin(), min(pl1$xMin(), pl2$xMin()))
  expect_equal(pl12$xMax(), max(pl1$xMax(), pl2$xMax()))
  expect_equal(pl12$xBy(), pl1$xBy())
  
  # discrete landscapes (require equal `xby`)
  pl1 <- pl_new(matrix(c(1, 3), ncol = 2L), exact = FALSE, xby = .2)
  pl2 <- pl_new(matrix(c(2, 4), ncol = 2L), exact = FALSE, xby = .2)
  pl12 <- pl1 + pl2
  expect_equal(pl12$xMin(), min(pl1$xMin(), pl2$xMin()))
  expect_equal(pl12$xMax(), max(pl1$xMax(), pl2$xMax()))
  expect_equal(pl12$xBy(), pl1$xBy())
  
  # exact and discrete landscapes (with incompatible `xby` in exact)
  pl1 <- pl_new(matrix(c(1, 3), ncol = 2L), exact = TRUE, xby = .1)
  pl2 <- pl_new(matrix(c(2, 4), ncol = 2L), exact = FALSE, xby = .2)
  pl12 <- pl1 + pl2
  expect_equal(pl12$xMin(), min(pl1$xMin(), pl2$xMin()))
  expect_equal(pl12$xMax(), max(pl1$xMax(), pl2$xMax()))
  expect_equal(pl12$xBy(), pl2$xBy())
  # TODO: Expose `epsi` to R and remove `round()`.
  expect_equal(unique(round(diff(pl12$getInternal()[, , 1L]), digits = 1L)), .2)
  expect_equal(
    pl12$getInternal()[, , 2L],
    # NOTE: assumes `xmin = 0`
    c(rep(0, 5L), seq(0, 1, .2), rep(1, 4), seq(1, 0, -.2))
  )
})
