
par(mfrow = c(2L, 2L), mar = rep(.5, 4))

# example PLs
x <- tdaunif::sample_klein_flat(60, ar = 2)
pd <- as_persistence(ripserr::vietoris_rips(x, dim = 1L, threshold = 2))
pl_e <- landscape(pd, degree = 1L, exact = TRUE)
pl_d <- landscape(pd, degree = 1L, exact = FALSE, xmax = 2, by = 0.05)
plot(pl_e, xaxt = "n", yaxt = "n")
plot(pl_d, xaxt = "n", yaxt = "n")

# vectorize PLs
vec_e <- pl_vectorize(pl_e)
vec_d <- pl_vectorize(pl_d)
length(vec_e)
attributes(vec_e)
length(vec_d)
attributes(vec_d)

# de-vectorized PL
pl_e_ <- pl_devectorize(vec_e)
plot(pl_e_, xaxt = "n", yaxt = "n")
pl_d_ <- pl_devectorize(vec_d)
plot(pl_d_, xaxt = "n", yaxt = "n")

par(mfrow = c(1L, 1L), mar = c(5.1, 4.1, 4.1, 2.1))