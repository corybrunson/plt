
set.seed(1)
pl_list_1 <- c()
for (i in seq(6)) {
 pc <- tdaunif::sample_lemniscate_gerono(n = 30, sd = 0.05)
 pd <- TDA::alphaComplexDiag(pc, maxdimension = 1)$diagram
 pd[,c(2,3)] <- sqrt(pd[,c(2,3)])
 pl <- pl_new(pd, degree = 1, xby = .01)
 pl_list_1 <- c(pl_list_1, pl)
}

set.seed(2)
pl_list_2 <- c()
for (i in seq(6)) {
 pc <- tdaunif::sample_circle(n = 30, sd = 0.05)
 pd <- TDA::alphaComplexDiag(pc, maxdimension = 1)$diagram
 pd[,c(2,3)] <- sqrt(pd[,c(2,3)])
 pl <- pl_new(pd, degree = 1, xby = .01)
 pl_list_2 <- c(pl_list_2, pl)
}

test_that("two-sided symmetry of pl_z_test() & pl_perm_test()", {
#pl_z_test
expect_equal(pl_z_test(pl_list_1,pl_list_2)$p.value, pl_z_test(pl_list_2,pl_list_1)$p.value)
#pl_perm_test
expect_equal(pl_perm_test(pl_list_1,pl_list_2)$p.value, pl_perm_test(pl_list_2,pl_list_1)$p.value, tolerance = 0.01)
})

 

# test one-sided symmetry of pl_z_test() & pl_perm_test()
