library(logPoolR)
source("R/nnn/config.R")
ms <- c(mu_1, mu_2, mu_3)
vs <- c(sigma_1, sigma_2, sigma_3) ^ 2
locking.weights <- c(w1, w2, w3)
thetastar <- logPoolR::pool_par_gauss(alpha = locking.weights,
                                      m = ms,
                                      v = vs)
###### Hyperparameter prep

mstar <- thetastar[1]
vstar <- thetastar[2] ^ 2

# mu_phi <- sin(mstar) * exp(-vstar/2) ## if functional(x) = sin(x)
mu_phi <- mstar ## functional(x) = x

k <- 2
mk <- ms[k]
vk <- vs[k]

ak <- (2 * vk - vstar) / (2 * vstar * vk)
bk <- (2 * vk * mstar - vstar * mk) / (2 * vstar * vk)
ck <- (2 * vk * mstar ^ 2 - vstar * mk ^ 2) / (2 * vstar * vk)

mm <- bk/ak
vv <- 1/(2*ak)

########

IS_kernel <- function(x)  functional(x) ^ 2 * exp(2 * dnorm(x, mstar, sqrt(vstar), log = TRUE) - dnorm(x, mk, sqrt(vk), log = TRUE))
var.quad <- integrate(IS_kernel, -Inf, Inf)$value - mu_phi ^ 2
var.quad

E_phi_sq <-  (vv + mm^2)
# E_phi_sq <- (exp(-2 * vv) * (exp(2 * vv) - cos(2 * mm))) / 2 ## if functional(x) = sin(x)

var.analytic <-
  exp(-ck + bk ^ 2 / ak)  * sqrt(vk/(2*ak)) / vstar * E_phi_sq - mu_phi ^ 2

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
mu_phi
var(I.hats) * M
var.quad
var.analytic
