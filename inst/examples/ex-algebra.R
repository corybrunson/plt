
par(mfrow = c(2L, 2L), mar = rep(.5, 4))

# example PLs
x <- tdaunif::sample_klein_flat(60, ar = 2)
pd <- ripserr::vietoris_rips(x, dim = 1L, threshold = 2)
pl_e <- pl_new(pd, degree = 1L, exact = TRUE)
pl_d <- pl_new(pd, degree = 1L, exact = FALSE, xmax = 2, xby = 0.05)
plot(pl_e, xaxt = "n", yaxt = "n")
plot(pl_d, xaxt = "n", yaxt = "n")

# vectorize PLs
vec_e <- pl_to_matrix(pl_e)
vec_d <- pl_to_matrix(pl_d)
length(vec_e)
attributes(vec_e)
length(vec_d)
attributes(vec_d)

# de-vectorized PL
pl_e_ <- pl_from_matrix(vec_e)
plot(pl_e_, xaxt = "n", yaxt = "n")
pl_d_ <- pl_from_matrix(vec_d)
plot(pl_d_, xaxt = "n", yaxt = "n")

par(mfrow = c(1L, 1L), mar = c(5.1, 4.1, 4.1, 2.1))

# vectorize a list
pl_lst <- list(pl_e, pl_discretize(pl_delimit(pl_e, xby = 0.05)), pl_d)
m <- pl_to_matrix(pl_lst)
# de-vectorize the matrix
pl_from_matrix(m)
pl_from_matrix(m, drop_levels = TRUE)
