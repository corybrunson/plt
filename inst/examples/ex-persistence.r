
# sample from a cylinder
x <- tdaunif::sample_cylinder_elliptical(n = 60)
# calculate PH using {ripserr}
ph <- ripserr::vietoris_rips(x, dim = 1)
as_persistence(ph)
# calculate PH using {TDA}
ph <- TDA::alphaComplexDiag(x, maxdimension = 1)
as_persistence(ph)
# create a 'persistence' object directly from a matrix
as_persistence(unclass(ph$diagram))
