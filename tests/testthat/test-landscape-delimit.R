
sq <- rbind(c(0, 0), c(0, 1), c(1, 1), c(1, 0))
ph <- TDA::alphaComplexDiag(sq)
# undo bug in `alphaComplexDiag()`
ph$diagram[, 2:3] <- sqrt(ph$diagram[, 2:3])
pl <- landscape(ph, degree = 1, xmax = 1, xby = 0.01)

test_that("`$delimit()` throws error iff limits contain not PL support", {
  supp <- pl_support(pl)
  # change `xmin`
  expect_no_error(pl_delimit(pl, xmin = supp[1L] - 0.01))
  expect_no_error(pl_delimit(pl, xmin = supp[1L]))
  expect_error(pl_delimit(pl, xmin = supp[1L] + 0.01), regexp = "support")
  # change `xmax`
  expect_no_error(pl_delimit(pl, xmax = supp[2L] + 0.01))
  expect_no_error(pl_delimit(pl, xmax = supp[2L]))
  expect_error(pl_delimit(pl, xmax = supp[2L] - 0.01), regexp = "support")
})

test_that("`$delimit()` preserves regularity", {
  lims <- pl_limits(pl)
  # change `xmin`
  pl1 <- pl_delimit(pl, xmin = lims[1L] - 0.01)
  expect_true(almostUnique(diff(pl1$getInternal()[, , 1L])))
  # change `xmax`
  pl2 <- pl_delimit(pl, xmax = lims[2L] + 0.01)
  expect_true(almostUnique(diff(pl2$getInternal()[, , 1L])))
})
