# This R function contains the code for computing the difference in RMST between the treatment and control groups
# for BayesDTE

library(splines2)
library(dplyr)

# Compare to integrate in R
surv_ctrl <- function(t, knots, bdry_knots, coef) {
  design <- iSpline(t, knots = knots, Boundary.knots = bdry_knots)
  exp(-design %*% coef)
}
rmst_ctrl_r <- function(t_upper, knots, bdry_knots, coef, lower = 0) {
  integrate(
    function(u) surv_ctrl(u, knots, bdry_knots, coef),
    lower = lower,
    upper = t_upper
  )$value
}
rmst_ctrl_r <- Vectorize(rmst_ctrl_r, "t_upper")

## Baseline hazard
haz_ctrl_r <- function(t, knots, bdry_knots, coef) {
  design <- iSpline(t, knots=knots, Boundary.knots = bdry_knots, derivs = 1)
  design %*% coef
}

## Hazard ratio between experimental arm and control arm
hazard_ratio_r <- function(t, tau, epsilon, hr_post) {
  g <- function(t) {
    alpha <- c((epsilon^3 - 3 * epsilon * tau^2 + hr_post * tau^2 * (3 * epsilon + 2 * tau) - 2 * tau^3) / epsilon^3,
               (6 * tau * (epsilon - hr_post * (epsilon + tau) + tau)) / epsilon^3,
               (3 * (-epsilon + hr_post * (epsilon + 2 * tau) - 2 * tau)) / epsilon^3,
               (2 * (1 - hr_post)) / epsilon^3)
    t_mat <- matrix(c(rep(1, length(t)), t, t^2, t^3), ncol = 4)
    t_mat %*% alpha
  }
  
  ratios <- c()
  ratios[which(t <= tau)]  <- 1
  ratios[which(t > tau & t <= tau + epsilon)]  <- g(t[which(t > tau & t <= tau + epsilon)])
  ratios[which(t > tau + epsilon)]  <- hr_post
  
  ratios
}

## Treatment group hazard
haz_trt_r <- function(t, knots, bdry_knots, coef, tau, epsilon, hr_post) {
  haz_ctrl_r(t, knots, bdry_knots, coef) * hazard_ratio_r(t, tau, epsilon, hr_post)
}

## Treatment group cumulative hazard: integrate numerically for transition window
cumhaz_trt_r <- function(t, knots, bdry_knots, coef, tau, epsilon, hr_post, lower = 0) {
  integrate(
    function(u) haz_trt_r(u, knots, bdry_knots, coef, tau, epsilon, hr_post),
    lower = lower,
    upper = t
  )$value
}
cumhaz_trt_r <- Vectorize(cumhaz_trt_r, "t")

surv_trt_r <- function(t, knots, bdry_knots, coef, tau, epsilon, hr_post, lower = 0) {
  exp(-cumhaz_trt_r(t, knots, bdry_knots, coef, tau, epsilon, hr_post))
}

rmst_trt_r <- function(t_upper, knots, bdry_knots, coef, tau, epsilon, hr_post, lower = 0) {
  integrate(
    function(u) surv_trt_r(u, knots, bdry_knots, coef, tau, epsilon, hr_post),
    lower = lower,
    upper = t_upper
  )$value
}
rmst_trt_r <- Vectorize(rmst_trt_r, "t_upper")