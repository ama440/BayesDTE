# This file contains the function to generate data from the DTE model using an
# i-spline on the cumulative hazard

# Load required libraries
library(splines2)

# Data generation
generate_data <- function(n, knots, bndry_knots, coef, treatment = FALSE, tau = NULL, 
                          epsilon = NULL, hr_post = NULL, seed = NULL) {
  ## Checks
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (treatment == TRUE & is.null(tau)) {
    stop("tau must be provided if treatment == TRUE")
  }
  
  if (treatment == TRUE & is.null(epsilon)) {
    stop("epsilon must be provided if treatment == TRUE")
  }
  
  if (treatment == TRUE & is.null(hr_post)) {
    stop("hr_post must be provided if treatment == TRUE")
  }
  
  # Check if number of knots matches number of coefficients
  if ((length(knots) + 4) != length(coef)) {
    stop("The length of coef must equal the number of knots.")
  }
  
  ## Baseline hazard
  haz_ctrl <- function(t, knots, bdry_knots, coef) {
    design <- iSpline(t, knots=knots, Boundary.knots = bdry_knots, derivs = 1)
    design %*% coef
  }
  
  ## Baseline cumulative hazard
  cumhaz_ctrl <- function(t, knots, bdry_knots, coef) {
    design <- iSpline(t, knots=knots, Boundary.knots = bdry_knots)
    design %*% coef
  }
  
  ## Baseline CDF
  cdf_ctrl <- function(t, knots, bdry_knots, coef) {
    1 - exp(-cumhaz_ctrl(t, knots, bdry_knots, coef))
  }
  
  ## Baseline Inverse CDF
  inverseCDF_ctrl <- function(u, knots, bdry_knots, coef) {
    uniroot(function(x) cdf_ctrl(x, knots, bdry_knots, coef) - u, 
            lower = 0, 
            upper = 1000)$root
  }
  inverseCDF_ctrl <- Vectorize(inverseCDF_ctrl, "u")
  
  ## Hazard ratio between experimental arm and control arm
  hazard_ratio <- function(t, tau, epsilon, hr_post) {
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
  haz_trt <- function(t, knots, bdry_knots, coef, tau, epsilon, hr_post) {
    haz_ctrl(t, knots, bdry_knots, coef) * hazard_ratio(t, tau, epsilon, hr_post)
  }
  
  ## Treatment group cumulative hazard: integrate numerically for transition window
  cumhaz_trt <- function(t, knots, bdry_knots, coef, tau, epsilon, hr_post, lower = 0) {
    integrate(
      function(u) haz_trt(u, knots, bdry_knots, coef, tau, epsilon, hr_post),
      lower = lower,
      upper = t
    )$value
  }
  cumhaz_trt <- Vectorize(cumhaz_trt, "t")
  
  
  ## Treatment group CDF
  cdf_trt <- function(t, knots, bdry_knots, coef, tau, epsilon, hr_post, lower = 0) {
    1 - exp(-cumhaz_trt(t, knots, bdry_knots, coef, tau, epsilon, hr_post, lower))
  }
  
  ## Treatment group Inverse CDF
  inverseCDF_trt <- function(u, knots, bdry_knots, coef, tau, epsilon, hr_post) {
    uniroot(function(x) cdf_trt(x, knots, bdry_knots, coef, tau, epsilon, hr_post) - u, 
            # lower = 1e-10, 
            lower = 0,
            upper = 1000, extendInt = "upX")$root
  }
  inverseCDF_trt <- Vectorize(inverseCDF_trt, "u")
  
  # Generate n uniform(0,1) samples
  u <- runif(n)
  
  # Generate from control group if treatment is false
  if (treatment == FALSE) {
    u[which(u > 0.99)] <- 0.99   # Make sure u isn't too large so that uniroot can find the value; these cases will be censored anyway so doesn't matter
    return(inverseCDF_ctrl(u, knots, bdry_knots, coef))
  }
  
  # Generate from treatment group if treatment is true
  if (treatment == TRUE) {
    u[which(u > 0.99)] <- 0.99   # Make sure u isn't too large so that uniroot can find the value; these cases will be censored anyway so doesn't matter
    return(inverseCDF_trt(u, knots, bdry_knots, coef, tau, epsilon, hr_post))
  }
}