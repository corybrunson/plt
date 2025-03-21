
# pl1s <- replicate(6, {
#   pc <- tdaunif::sample_lemniscate_gerono(60, sd = .1)
#   pd <- ripserr::vietoris_rips(pc, dim = 1L, threshold = 2, p = 2)
#   pl_new(pd, degree = 1, xby = .01)
# })
# # samples of landscapes from circles
# pl2s <- replicate(8, {
#   pc <- tdaunif::sample_circle(60, sd = .1) / 2
#   pd <- ripserr::vietoris_rips(pc, dim = 1L, threshold = 2, p = 2)
#   pl_new(pd, degree = 1, xby = .01)
# })
# 
# 
# test_that("pd_z_test() successfully runs",{
#   expect_no_error(pd_z_test(pl1s,pl2s))
# })
# 






###############################################################################

# set.seed(1)
# pl_list_1 <- c()
# for (i in seq(6)) {
#  pc <- tdaunif::sample_lemniscate_gerono(n = 60, sd = 0.05)
#  pd <- TDA::alphaComplexDiag(pc, maxdimension = 1)$diagram
#  pd[,c(2,3)] <- sqrt(pd[,c(2,3)])
#  pl <- pl_new(pd, degree = 1, xby = .01)
#  pl_list_1 <- c(pl_list_1, pl)
# }
# set.seed(2)
# pl_list_2 <- c()
# for (i in seq(6)) {
#  pc <- tdaunif::sample_circle(n = 60, sd = 0.05)
#  pd <- TDA::alphaComplexDiag(pc, maxdimension = 1)$diagram
#  pd[,c(2,3)] <- sqrt(pd[,c(2,3)])
#  pl <- pl_new(pd, degree = 1, xby = .01)
#  pl_list_2 <- c(pl_list_2, pl)
# }

#test_that("symmetry of 3 statistical functions", {
  #pd_z_test
  #expect_equal(pd_z_test(pl_list_1,pl_list_2),pd_z_test(pl_list_2,pl_list_1))
  #pl_z_test
  #expect_equal(pl_z_test(pl_list_1,pl_list_2),pl_z_test(pl_list_2,pl_list_1))
  #pl_perm_test
  #expect_equal(pl_perm_test(pl_list_1,pl_list_2),pl_perm_test(pl_list_2,pl_list_1))
#})
