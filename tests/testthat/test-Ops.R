
test_that("Addition expands domain as needed.", {
  
  # exact landscapes
  pl1 <- landscape(matrix(c(1, 3), ncol = 2L), exact = TRUE, by = .1)
  pl2 <- landscape(matrix(c(2, 4), ncol = 2L), exact = TRUE, by = .01)
  pl12 <- pl1 + pl2
  expect_equal(pl12$xMin(), min(pl1$xMin(), pl2$xMin()))
  expect_equal(pl12$xMax(), max(pl1$xMax(), pl2$xMax()))
  expect_equal(pl12$xBy(), pl1$xBy())
  
  # discrete landscapes (require equal `by`)
  pl1 <- landscape(matrix(c(1, 3), ncol = 2L), exact = FALSE, by = .1)
  pl2 <- landscape(matrix(c(2, 4), ncol = 2L), exact = FALSE, by = .1)
  pl12 <- pl1 + pl2
  expect_equal(pl12$xMin(), min(pl1$xMin(), pl2$xMin()))
  expect_equal(pl12$xMax(), max(pl1$xMax(), pl2$xMax()))
  expect_equal(pl12$xBy(), pl1$xBy())
  
  # exact and discrete landscapes
  pl1 <- landscape(matrix(c(1, 3), ncol = 2L), exact = TRUE, by = .1)
  pl2 <- landscape(matrix(c(2, 4), ncol = 2L), exact = FALSE, by = .2)
  pl12 <- pl1 + pl2
  expect_equal(pl12$xMin(), min(pl1$xMin(), pl2$xMin()))
  expect_equal(pl12$xMax(), max(pl1$xMax(), pl2$xMax()))
  expect_equal(pl12$xBy(), pl2$xBy())
  # TODO: Expose `epsi` to R and remove `round()`.
  expect_equal(unique(round(diff(pl12$getInternal()[, , 1L]), digits = 1L)), .2)
})
