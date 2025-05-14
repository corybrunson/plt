
# multiply by an indicator function
pd <- data.frame(dim = 1, birth = c(0, 1), death = c(2, 4))
f <- list(c(.5, 3), c(1.5, 2), c(1.75, 2))

# exact landscape
pl_e <- pl_new(pd, degree = 1, exact = TRUE)
par(mfcol = c(2L, 2L), mar = c(0, 2, 0, 2))
plot(pl_e, xlim = c(0, 5))
# pl_ef <- pl$indicator(f, 0)
pl_ef <- pl_indicator(pl_e, f)
plot(pl_ef, xlim = c(0, 5))
# discrete landscape
pl_d <- pl_new(pd, degree = 1, xmin = 0, xmax = 5, xby = .2)
plot(pl_d)
# pl_df <- pl_d$indicator(f, 0)
pl_df <- pl_indicator(pl_d, f)
plot(pl_df, xlim = c(0, 5))
par(mfrow = c(1L, 1L), mar = c(5.1, 4.1, 4.1, 2.1))

# integrate PL multiplied by an indicator function
pl_indicator_form(pl_d, supports = f, p = 5)
pl_integrate(pl_df, p = 5)
pl_indicator_form(pl_e, supports = f, p = 5)
pl_integrate(pl_ef, p = 5)
