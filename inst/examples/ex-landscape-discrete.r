
# sample points
points <- tdaunif::sample_torus_tube(60L, 2.5)

# compute persistence data, retaining parameters
pd <- as_persistence(TDA::ripsDiag(points, maxdimension = 2L, maxscale = 3))
print(dim(pd$pairs[[2L]]))

# compute persistence landscape for 1-dimensional cycles
pl <- landscape(pd, degree = 1L)
print(pl)

# landscape dimensions
print(dim(pl$getInternal()))
# landscape values
print(head(pl$getInternal()))
# plot landscape
plot(pl)

# custom parameters
pl <- landscape(pd, degree = 1L, by = 0.1, xmax = 2)
print(pl)
plot(pl)
