
# sample points
points <- tdaunif::sample_torus_tube(100, 5)

# compute persistent homology
(pd <- ripserr::vietoris_rips(points, dim = 2L, threshold = 1))

# compute persistence landscapes for 0-cycles
(pl <- pl_new(pd, degree = 1, exact = TRUE))

# first landscape layer
print(pl$getInternal()[[1L]])
# plot all landscape layers
plot(pl)

# coerce to discrete at different resolutions
plot(pl, xlim = c(0, .5))
pl <- pl_delimit(pl, xmin = 0, xmax = .5, xby = .01)
pl_ <- pl_discretize(pl)
plot(pl_)
pl <- pl_delimit(pl, xmin = 0, xmax = .5, xby = .001)
pl_ <- pl_discretize(pl)
plot(pl_)

# ensure grid when discretizing
pl_cut <- pl_discretize(pl_delimit(pl, xby = 0.1))
pl_cut$getInternal()[1, , 1]
pl_cut <- pl_discretize(pl_delimit(pl, xmin = 0, xmax = pi, xby = 0.1))
pl_cut$getInternal()[1, , 1]
