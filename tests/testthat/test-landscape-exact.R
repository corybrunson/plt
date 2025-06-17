test_that("PL is correct for one persistence pair.", {
  p1 <- matrix(c(0, 2), nrow = 1, ncol = 2)
  pd <- p1
  pl <- pl_new(pd, exact = TRUE)
  expected <- list(matrix(c(-Inf, 0, 1, 2, Inf, 0, 0, 1, 0, 0),
                          nrow = 5, ncol = 2))
  
  expect_equal(pl$getInternal(), expected, ignore_attr = FALSE)
})

scale <- function(h, pl_exact){
  out <- list()
  
  for(i in 1:length(pl_exact)){
    level <- pl_exact[[i]]
    x <- level[,1]
    y <- h*level[,2]
    out <- list(out,c(x,y))
  }
  
  return(out)
}

test_that("PL sum is correct for simple case.", {
  pd1 <- matrix(c(0, 2), nrow = 1, ncol = 2)
  pd2 <- matrix(c(1, 2), nrow = 1, ncol = 2)
  
  pl1 <- pl_new(pd1, exact = TRUE)
  pl2 <- pl_new(pd2, exact = TRUE)
  
  pl <- pl_sum(list(pl1,pl2))
  expected <- list(matrix(c(-Inf, 0, 1, 1.5, 2, Inf, 0, 0, 1, 1, 0, 0),
                          nrow = 6, ncol = 2))
  
  expect_equal(pl$getInternal(), expected, ignore_attr = FALSE)
})

test_that("add PL is correct for simple case.", {
  pd1 <- matrix(c(0, 2), nrow = 1, ncol = 2)
  pd2 <- matrix(c(1, 2), nrow = 1, ncol = 2)
  
  pl1 <- pl_new(pd1, exact = TRUE)
  pl2 <- pl_new(pd2, exact = TRUE)
  
  pl <- pl1$add(pl2)
  expected <- list(matrix(c(-Inf, 0, 1, 1.5, 2, Inf, 0, 0, 1, 1, 0, 0),
                          nrow = 6, ncol = 2))
  
  expect_equal(pl$getInternal(), expected, ignore_attr = FALSE)
})

x <- tdaunif::sample_circle(100)
pd <- as_persistence(ripserr::vietoris_rips(x, max_dim = 1L, threshold = 2))

test_that("`discretize` from exact is correct", {
  pl <- pl_new(pd$pairs[[1]], exact=TRUE)
  
  expect_error(pl$discretize(), NA)
})

test_that("getInternal from discrete is correct from diagram", {
  pd <- suppressWarnings(
    as_persistence(ripserr::vietoris_rips(x, max_dim = 1L, threshold = 2)))
  pl <- pl_new(pd, degree = 1L, exact = TRUE,
                  xmax = 2.5, xby = 0.1)
  
  pdref <- suppressWarnings(
    as_persistence(ripserr::vietoris_rips(x, max_dim = 1L, threshold = 2)))
  plref <- pl_new(pdref, degree = 1L, exact = TRUE,
                  xmax = 2.5, xby = 0.1)
  
  expect_equal(pl$getInternal(), plref$getInternal())
})
