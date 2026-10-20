## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(fig.width = 6, fig.height = 4, collapse = TRUE, comment = "#>")
set.seed(20260925)

## -----------------------------------------------------------------------------
library(qdfauction)
data(timber)
str(timber)

## -----------------------------------------------------------------------------
b3  <- as.numeric(timber$n3)
# bandwidth = "bcv" selects h by biased cross-validation (about 20 s here);
# we fix h at the value it selects to keep the vignette fast
q_ip <- qdf(b3, method = "indirect_poisson", bandwidth = 0.02139)
q_b  <- qdf(b3, method = "bernstein",        bandwidth = "bcv")
q_w  <- qdf(b3, method = "wavelet_block",    bandwidth = 0.15)
q_ip
plot(q_ip, ylim = c(0, 400), main = "Quantile density of bids, N = 3")
lines(q_b$u, q_b$q, col = 2); lines(q_w$u, q_w$q, col = 4)
legend("topleft", c("indirect Poisson", "Bernstein", "block wavelet"),
       col = c(1, 2, 4), lty = 1, bty = "n")

## -----------------------------------------------------------------------------
fv_ip  <- fpa_values(timber$n3, method = "indirect_poisson", bandwidth = q_ip$h,
                     monotone = "spline")           # reuse the BCV bandwidth from above
fv_gpv <- fpa_values(timber$n3, method = "gpv")
fv_ip
plot(fv_ip, main = "Inverse bid function, N = 3")
points(fv_gpv$bids, fv_gpv$values, col = 2, pch = 3, cex = 0.6)
legend("topleft", c("indirect Poisson (raw)", "monotone fit", "GPV"),
       col = c("grey50", 4, 2), pch = c(20, NA, 3), lty = c(NA, 1, NA), bty = "n")

## -----------------------------------------------------------------------------
fv2 <- fpa_values(timber$n2, method = "indirect_poisson", bandwidth = 0.02108,   # BCV choice
                  monotone = "spline")
rbind(`N = 2` = quantile(fv2$values - fv2$bids, c(.25, .5, .75)),
      `N = 3` = quantile(fv_ip$values - fv_ip$bids, c(.25, .5, .75)))

## -----------------------------------------------------------------------------
fv_ra <- fpa_values(timber$n3, method = "indirect_poisson", bandwidth = q_ip$h,
                    eta = 0.3, monotone = "spline")
median(fv_ra$values - fv_ra$bids) / median(fv_ip$values - fv_ip$bids)   # = 1 - eta

## -----------------------------------------------------------------------------
fd_w  <- fpa_density(timber$n3, method = "wavelet_ll", j0 = 4, h = 0.15)
fd_hh <- fpa_density(timber$n3, method = "kde_hh")
fd_ms <- fpa_density(timber$n3, method = "marmer_shneyerov")
fd_q  <- fpa_density(timber$n3, method = "wavelet_ll", first_stage = fv_ip, j0 = 4, h = 0.15)
fd_w
plot(fd_w, xlim = c(0, 400), ylim = c(0, 0.02), main = "Private-value density, N = 3")
lines(fd_hh$x, fd_hh$f, col = 2, lty = 2)
lines(fd_ms$x, fd_ms$f, col = 3, lty = 3)
lines(fd_q$x,  fd_q$f,  col = 4, lty = 4)
legend("topright", c("wavelet LL (GPV-HH first stage)", "GPV-HH kernel",
                     "Marmer-Shneyerov", "wavelet LL (QDF first stage)"),
       col = 1:4, lty = 1:4, bty = "n")

## ----eval = FALSE-------------------------------------------------------------
# fd_cv <- fpa_density(timber$n3, method = "wavelet_ll", j0 = "cv")
# fpa_density_cv(timber$n3, methods = c("wavelet_ll", "kde_hh", "marmer_shneyerov"),
#                folds = 5, j0 = fd_cv$j0, h = fd_cv$h)

## -----------------------------------------------------------------------------
ci3 <- fpa_density_ci(timber$n3, B = 199)
ci3
plot(ci3, main = "GPV density with 95% uniform band, N = 3")

## -----------------------------------------------------------------------------
r  <- seq(20, 120, by = 5)
R3 <- fpa_revenue(fv_ip, r, v_s = 20)
opt <- fpa_optimal_reserve(fv_ip, v_s = 20)
plot(r, R3, type = "l", lwd = 2, xlab = "reserve price ($/unit)",
     ylab = "expected revenue", main = "Expected revenue, N = 3")
abline(v = opt$r_star, lty = 2)
opt[c("r_star", "revenue")]

## -----------------------------------------------------------------------------
rb <- fpa_revenue_boot(timber$n3, r = seq(30, 100, by = 10), B = 30, v_s = 20,
                       method = "bernstein", bandwidth = "wbcv", m_range = c(10, 200))
rb
plot(rb, main = "Expected revenue with 95% pointwise CI and uniform band")

## -----------------------------------------------------------------------------
pdf <- function(x) dgamma(x, 10, 2); cdf <- function(x) pgamma(x, 10, 2)
mc <- fpa_montecarlo(R = 3, L = 20, n = 5, Q = function(u) qgamma(u, 10, 2), seed = 1,
  estimators = list(
    B_B  = function(b, v) fpa_values(b, method = "bernstein",        bandwidth = "bcv"),
    B_IP = function(b, v) fpa_values(b, method = "indirect_poisson", bandwidth = "bcv"),
    Conv = function(b, v) fpa_values(b, method = "gpv")),
  score = function(fit, b, v) fpa_evaluate_values(fit, pdf, cdf, v))
aggregate(cbind(mise, loglik, ks) ~ estimator, mc, mean)

