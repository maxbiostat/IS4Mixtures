library(logPoolR)
source("R/nnn/config.R")
ms <- c(mu_1, mu_2, mu_3)
vs <- c(sigma_1, sigma_2, sigma_3) ^ 2
locking.weights <- c(w1, w2, w3)
thetastar <- logPoolR::pool_par_gauss(alpha = locking.weights,
                                      m = ms,
                                      v = vs)
mstar <- thetastar[1]
vstar <- thetastar[2] ^ 2
k <- 1
mk <- ms[k]
vk <- vs[k]

ak <- (2 * vk - vstar) / (2 * vstar * vk)
bk <- (2 * vk * mstar - vstar * mk) / (2 * vstar * vk)
ck <- (2 * vk * mstar ^ 2 - vstar * mk ^ 2) / (2 * vstar * vk)


IS_kernel <- function(x)  functional(x) ^ 2 * exp(2 * dnorm(x, mstar, sqrt(vstar), log = TRUE) - dnorm(x, mk, sqrt(vk), log = TRUE))
var.quad <- integrate(IS_kernel, -Inf, Inf)$value - mstar ^ 2
var.quad

var.analytic <-
  exp(-ck + bk ^ 2 / ak)  * sqrt(vk/(2*ak)) / vstar * (1 / (2*ak) + (bk / ak)^2) - mstar ^ 2

get_est <- function(xs) {
  vals <-
    functional(xs) * exp(dnorm(xs, mstar, sqrt(vstar), log = TRUE) - dnorm(xs, mk, sqrt(vk), log = TRUE))
  return(mean(vals))
}

M <- 1E4
Nrep <- 500
simus <- matrix(rnorm(M * Nrep, mk, sqrt(vk)),
                ncol = M, nrow = Nrep)

I.hats <- apply(simus, 1, get_est)

hist(I.hats)
mean(I.hats)
mstar
var(I.hats) * M
var.quad
var.analytic
