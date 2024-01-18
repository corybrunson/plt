
# sample points
points <- tdaunif::sample_torus_tube(100, 5)

# compute persistence
pd <- as_persistence(ripserr::vietoris_rips(points, dim = 2L, threshold = 1))
print(pd)

# compute persistence landscapes for 0-cycles
pl <- landscape(pd, degree = 1, exact = TRUE)
print(pl)

# first landscape layer
print(pl$getInternal()[[1L]])
# plot all landscape layers
plot(pl)

# # coerce to discrete
# pl_ <- pl_discretize(pl, xmin = 0, xmax = .5, by = .005)
# plot(pl_)
