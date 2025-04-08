
test_that("limits that don't contain support generate a warning/error", {
  sq <- rbind(c(0, 0), c(0, 1), c(1, 1), c(1, 0))
  ph <- TDA::ripsDiag(sq, maxdimension = 2L, maxscale = 1.5)
  expect_error(pl_new(ph, degree = 1, xmax = .4))
  expect_error(pl_new(ph, degree = 1, xmin = .3, xmax = 1))
})

test_that("PL is correct for one persistence pair.", {
  p1 <- matrix(c(0, 2), nrow = 1L, ncol = 2L)
  pd <- p1
  pl <- pl_new(pd, exact = FALSE, xmin = 0, xmax = 5, xby = 0.1)
  
  x_val <- seq(0, 5, 0.1)
  y_val_0 <- seq(0, 1, 0.1)
  y_val_1 <- seq(1-0.1, 0, -0.1)
  y_zeros <- rep(0, length(x_val) - (length(y_val_0) + length(y_val_1)))
  y_val <- c(y_val_0, y_val_1, y_zeros)
  expected <- array(cbind(x_val, y_val), c(1L, 51L, 2L))
  
  expect_equal(pl$getInternal(), expected, ignore_attr = FALSE)
})

x <- tdaunif::sample_circle(100L)
y <- tdaunif::sample_circle(100L)
pd <- as_persistence(ripserr::vietoris_rips(x, max_dim = 1L, threshold = 2))
pd2 <- as_persistence(ripserr::vietoris_rips(y, max_dim = 1L, threshold = 2))

matchDimension <- function(pl_1, pl_2){
  if(dim(pl_1)[1] < dim(pl_2)[1]){
    pl_1_t <- array(0, dim(pl_2))
    pl_1_t[1:dim(pl_1)[1],,] = pl_1
    pl_1_t[,,1] = pl_2[,,1]
    pl_1 = pl_1_t
  }
  
  else{
    pl_2_t <- array(0, dim(pl_1))
    pl_2_t[1:dim(pl_2)[1],,] = pl_2
    pl_2_t[,,1] = pl_1[,,1]
    pl_2 = pl_2_t
  }
  
  list(pl_1,pl_2)
  
}

sumPL <- function(pl1,pl2){
  pl_l = matchDimension(pl1, pl2)
  pl_1 = pl_l[[1]]
  pl_2 = pl_l[[2]]
  
  x <- pl_1[,,1]
  y <- (pl_1[,,2]+pl_2[,,2])
  
  array(c(x,y), dim=dim(pl_1))
}

scalePL <- function(lambda, pl){
  x <- pl[,,1]
  y <- lambda*pl[,,2]
  
  array(c(x,y), dim=dim(pl))
}

innerPL <- function(pl1,pl2){
  pl_l <- matchDimension(pl1,pl2)
  pl_1 <- as.vector(pl_l[[1]][,,2])
  pl_2 <- as.vector(pl_l[[2]][,,2])
  
  
  inner <- sum(pl_1%*%pl_2)
  inner*(pl1[1,2,1]-pl1[1,1,1])
}

test_that("PL sum is correct.", {
  pl <- pl_new(pd$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
  pl2 <- pl_new(pd2$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
  
  pl_d = pl$getInternal()
  pl2_d = pl2$getInternal()
  
  expect_equal(pl_sum(list(pl, pl2))$getInternal(), sumPL(pl_d, pl2_d))
})


test_that("add PL is correct.", {
  pl <- pl_new(pd$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
  pl2 <- pl_new(pd2$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
  
  pl_d = pl$getInternal()
  pl2_d = pl2$getInternal()
  
  expect_equal(pl$add(pl2)$getInternal(), sumPL(pl_d, pl2_d))
})


# test_that("PL scale is correct.", {
#   pl <- pl_new(pd$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
#   pl_d = pl$getInternal()
#   
#   expect_equal(PLscale(0.5, pl)$getInternal(), scalePL(0.5, pl_d))
# })

test_that("scale PL is correct.", {
  pl <- pl_new(pd$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
  pl_d = pl$getInternal()
  
  expect_equal(pl$scale(0.5)$getInternal(), scalePL(0.5, pl_d))
})

test_that("average PL is correct.", {
  pl <- pl_new(pd$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
  pl2 <- pl_new(pd2$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
  
  pl_d = pl$getInternal()
  pl2_d = pl2$getInternal()
  
  expect_equal(scalePL(0.5, sumPL(pl_d, pl2_d)),
               pl_mean(list(pl,pl2))$getInternal())
})

# test_that("PL inner is correct.", {
#   pl <- pl_new(pd$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
#   pl2 <- pl_new(pd2$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
#   
#   pl_d = pl$getInternal()
#   pl2_d = pl2$getInternal()
#   
#   expect_equal(PLinner(pl,pl2), innerPL(pl_d,pl2_d))
# })

test_that("inner PL is correct.", {
  pl <- pl_new(pd$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
  pl2 <- pl_new(pd2$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
  
  pl_d = pl$getInternal()
  pl2_d = pl2$getInternal()
  
  expect_equal(pl$inner(pl2), innerPL(pl_d,pl2_d))
})

test_that("toExact from discrete is correct", {
  pl <- pl_new(pd$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
  
  expect_error(pl$toExact())
})

test_that("`discretize` from discrete is correct", {
  pl <- pl_new(pd$pairs[[1]], exact=FALSE, xmax=2.5, xby=0.1)
  
  expect_warning(pl$discretize(), "discrete")
})

test_that("getInternal from discrete is correct from diagram", {
  pd <- suppressWarnings(
    as_persistence(ripserr::vietoris_rips(x, dim = 1L, threshold = 2)))
  pl <- pl_new(pd, degree = 1, exact = FALSE,
                  xmax = 2.5, xby = 0.1)
  
  pdref <- suppressWarnings(
    as_persistence(ripserr::vietoris_rips(x, dim = 1L, threshold = 2)))
  plref <- pl_new(pdref, degree = 1, exact = FALSE,
                     xmax = 2.5, xby = 0.1)
  
  expect_equal(pl$getInternal(), plref$getInternal())
})

