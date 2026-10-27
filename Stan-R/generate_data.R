# This file contains the function to generate data from the DTE model using an
# i-spline on the cumulative hazard

# Load required libraries
library(splines2)

# Data generation
generate_data <- function(n, knots, bdry_knots, coef, treatment = FALSE, tau = NULL, 
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

# ------------------------------------------------------------------------------
# Data generation under parametric baseline hazards
#
# The functions below generate data from the DTE model using a piecewise
# exponential (exponential as a special case), Weibull, or log-logistic baseline
# hazard for the control arm. The hazard ratio between the treatment and control
# arms is the same as in generate_data(): 1 up to tau, a cubic transition on
# (tau, tau + epsilon], and hr_post afterwards. Times are drawn by inverting the
# cumulative hazard at -log(U), using the closed-form inverse of the baseline
# cumulative hazard outside the transition window.
# ------------------------------------------------------------------------------

## Hazard ratio between experimental arm and control arm (same as in generate_data())
hazard_ratio_dte <- function(t, tau, epsilon, hr_post) {
  g <- function(t) {
    alpha <- c((epsilon^3 - 3 * epsilon * tau^2 + hr_post * tau^2 * (3 * epsilon + 2 * tau) - 2 * tau^3) / epsilon^3,
               (6 * tau * (epsilon - hr_post * (epsilon + tau) + tau)) / epsilon^3,
               (3 * (-epsilon + hr_post * (epsilon + 2 * tau) - 2 * tau)) / epsilon^3,
               (2 * (1 - hr_post)) / epsilon^3)
    t_mat <- matrix(c(rep(1, length(t)), t, t^2, t^3), ncol = 4)
    as.vector(t_mat %*% alpha)
  }
  
  ratios <- rep(NA_real_, length(t))
  ratios[which(t <= tau)]  <- 1
  ratios[which(t > tau & t <= tau + epsilon)]  <- g(t[which(t > tau & t <= tau + epsilon)])
  ratios[which(t > tau + epsilon)]  <- hr_post
  
  ratios
}

## Generic generator given the baseline hazard, cumulative hazard, and inverse
## cumulative hazard of the control arm
generate_data_haz <- function(n, haz_ctrl, cumhaz_ctrl, inv_cumhaz_ctrl, 
                              treatment = FALSE, tau = NULL, epsilon = NULL, 
                              hr_post = NULL, seed = NULL) {
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
  
  # Target cumulative hazard values: T = H^{-1}(-log(U))
  target <- -log(runif(n))
  
  # Generate from control group if treatment is false
  if (treatment == FALSE) {
    return(inv_cumhaz_ctrl(target))
  }
  
  # Generate from treatment group if treatment is true
  ## Treatment cumulative hazard accumulated over the transition window (tau, t]
  window_cumhaz <- function(t) {
    if (epsilon <= 0 || t <= tau) return(0)
    integrate(function(u) haz_ctrl(u) * hazard_ratio_dte(u, tau, epsilon, hr_post),
              lower = tau, upper = t)$value
  }
  
  H_tau <- cumhaz_ctrl(tau)                  # Treatment cumulative hazard at tau
  H_end <- H_tau + window_cumhaz(tau + epsilon) # Treatment cumulative hazard at tau + epsilon
  
  times <- numeric(n)
  
  # Before tau: same as control
  pre <- target <= H_tau
  times[pre] <- inv_cumhaz_ctrl(target[pre])
  
  # Transition window: solve numerically
  mid <- which(target > H_tau & target <= H_end)
  times[mid] <- vapply(target[mid], function(h) {
    uniroot(function(x) H_tau + window_cumhaz(x) - h,
            lower = tau, upper = tau + epsilon)$root
  }, numeric(1))
  
  # After tau + epsilon: proportional hazards with hr_post
  post <- target > H_end
  times[post] <- inv_cumhaz_ctrl(cumhaz_ctrl(tau + epsilon) + (target[post] - H_end) / hr_post)
  
  times
}

## Piecewise exponential baseline hazard
# rates: hazard on each interval; breaks: interior cut points (length(rates) - 1),
# so intervals are [0, breaks[1]), [breaks[1], breaks[2]), ..., [breaks[k], Inf).
# A single rate with breaks = NULL gives the exponential distribution.
generate_data_pwexp <- function(n, rates, breaks = NULL, treatment = FALSE, tau = NULL, 
                                epsilon = NULL, hr_post = NULL, seed = NULL) {
  if (any(rates <= 0)) {
    stop("rates must be positive")
  }
  
  if (length(breaks) != length(rates) - 1) {
    stop("The length of breaks must be one less than the length of rates.")
  }
  
  if (length(breaks) > 0 && (any(breaks <= 0) || is.unsorted(breaks, strictly = TRUE))) {
    stop("breaks must be positive and strictly increasing")
  }
  
  cuts <- c(0, breaks)
  cumhaz_cuts <- c(0, cumsum(rates[-length(rates)] * diff(cuts)))  # Cumulative hazard at each cut
  
  haz_ctrl <- function(t) {
    rates[findInterval(t, cuts)]
  }
  
  cumhaz_ctrl <- function(t) {
    j <- findInterval(t, cuts)
    cumhaz_cuts[j] + rates[j] * (t - cuts[j])
  }
  
  inv_cumhaz_ctrl <- function(h) {
    j <- findInterval(h, cumhaz_cuts)
    cuts[j] + (h - cumhaz_cuts[j]) / rates[j]
  }
  
  generate_data_haz(n, haz_ctrl, cumhaz_ctrl, inv_cumhaz_ctrl, treatment = treatment, 
                    tau = tau, epsilon = epsilon, hr_post = hr_post, seed = seed)
}

## Weibull baseline hazard
# Same parameterization as rweibull(): h(t) = (shape / scale) * (t / scale)^(shape - 1),
# H(t) = (t / scale)^shape
generate_data_weibull <- function(n, shape, scale, treatment = FALSE, tau = NULL, 
                                  epsilon = NULL, hr_post = NULL, seed = NULL) {
  if (shape <= 0 || scale <= 0) {
    stop("shape and scale must be positive")
  }
  
  haz_ctrl <- function(t) {
    (shape / scale) * (t / scale)^(shape - 1)
  }
  
  cumhaz_ctrl <- function(t) {
    (t / scale)^shape
  }
  
  inv_cumhaz_ctrl <- function(h) {
    scale * h^(1 / shape)
  }
  
  generate_data_haz(n, haz_ctrl, cumhaz_ctrl, inv_cumhaz_ctrl, treatment = treatment, 
                    tau = tau, epsilon = epsilon, hr_post = hr_post, seed = seed)
}

## Log-logistic baseline hazard
# S(t) = 1 / (1 + (t / scale)^shape), so scale is the median survival time,
# h(t) = (shape / scale) * (t / scale)^(shape - 1) / (1 + (t / scale)^shape),
# H(t) = log(1 + (t / scale)^shape)
generate_data_loglogistic <- function(n, shape, scale, treatment = FALSE, tau = NULL, 
                                      epsilon = NULL, hr_post = NULL, seed = NULL) {
  if (shape <= 0 || scale <= 0) {
    stop("shape and scale must be positive")
  }
  
  haz_ctrl <- function(t) {
    (shape / scale) * (t / scale)^(shape - 1) / (1 + (t / scale)^shape)
  }
  
  cumhaz_ctrl <- function(t) {
    log1p((t / scale)^shape)
  }
  
  inv_cumhaz_ctrl <- function(h) {
    scale * expm1(h)^(1 / shape)
  }
  
  generate_data_haz(n, haz_ctrl, cumhaz_ctrl, inv_cumhaz_ctrl, treatment = treatment, 
                    tau = tau, epsilon = epsilon, hr_post = hr_post, seed = seed)
}
