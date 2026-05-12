library(cmdstanr)
library(rstan)
library(posterior)
library(statmod)
library(splines2)

# Load cmdstan model
bayes.dte <- cmdstanr::cmdstan_model("Stan-R/dte_pp.stan")

# Load R functions for baseline and experimental arm cumulative hazards
source("Stan-R/rmst_funs.R")

#' Get stan data for BayesDTE
#' 
#' Obtains stan data for BayesDTE
#' 
#' @param t vector of observed failure times in current data
#' @param delta vector of event indicators in current data
#' @param z vector of treatment indicators in current data
#' @param t_hist vector of observed failure times in historical data (Optional)
#' @param delta_hist vector of event indicators in historical data (Optional)
#' @param z_hist vector of treatment indicators in historical data (Optional)
#' @param knots knot locations
#' @param n_nodes number of nodes to use for Gaussian quadrature
#' @param tau_val
#' @param epsilon_val
#' @param a0
#' @param prior_tau character: "uniform" (default), "mixture", or "half-normal" prior on tau
#' @param mixture_weight numeric: probability weight for the point mass spike at 0 (used if prior_tau == "mixture")
#' @param sigma_tau numeric: standard deviation for half-normal prior (used if prior_tau == "half-normal")
#'
#' @return a list giving data to be passed onto `rstan` or `cmdstanr` models.
get.standata.dte <- function(
    t, delta, z, 
    t_hist = numeric(0), delta_hist = numeric(0), z_hist = numeric(0), 
    knots, n_nodes, 
    tau_val = numeric(0), tau_max = numeric(0), epsilon_val = numeric(0), 
    sd_epsilon = 10, rho = 0, sd_beta = 10, 
    t_upper, a0 = 0,
    prior_tau = c("uniform", "mixture", "half-normal"), 
    mixture_weight = 0.5, sigma_tau = 1.0) {
  
  # Validate string input and ensure it matches one of the expected options
  prior_tau <- match.arg(prior_tau)
  
  # Map string to the integer required by Stan
  prior_tau_int <- switch(prior_tau,
                          "uniform" = 0,
                          "mixture" = 1,
                          "half-normal" = 2)
  
  n <- length(t)
  df <- length(knots) + 4
  t_max <- max(t)
  bdry_knots <- c(0, t_max)
  
  # Compute Gaussian quadrature
  gq <- gauss.quad(n = n_nodes, kind = "legendre")
  nodes <- gq$nodes
  weights <- gq$weights
  
  # Split current data
  t00 <- t[delta == 0 & z == 0]
  t01 <- t[delta == 0 & z == 1]
  t10 <- t[delta == 1 & z == 0]
  t11 <- t[delta == 1 & z == 1]
  
  # Safe spline constructor for when vectors are length 0
  safe_ispline <- function(x, derivs = 0) {
    if (length(x) == 0) return(matrix(0, nrow = 0, ncol = df))
    mat <- as.matrix(iSpline(x, knots = knots, Boundary.knots = bdry_knots, derivs = derivs))
    mat <- pmax(mat, 0)
    return(mat)
  }
  
  # Handle historical data (defaults to empty vectors if a0 == 0 or missing)
  if (a0 > 0 && length(t_hist) > 0) {
    t00_hist <- t_hist[delta_hist == 0 & z_hist == 0]
    t01_hist <- t_hist[delta_hist == 0 & z_hist == 1]
    t10_hist <- t_hist[delta_hist == 1 & z_hist == 0]
    t11_hist <- t_hist[delta_hist == 1 & z_hist == 1]
  } else {
    t00_hist <- numeric(0)
    t01_hist <- numeric(0)
    t10_hist <- numeric(0)
    t11_hist <- numeric(0)
  }
  n_hist <- length(t00_hist) + length(t01_hist) + length(t10_hist) + length(t11_hist)
  
  # Construct hazard design matrix for GQ points, censored and in treatment group
  des_gq_array_01 <- array(NA_real_, dim = c(n_nodes, length(t01), df))
  for (l in seq_len(n_nodes)) {
    input <- t01/2 * (nodes[l] + 1)
    des_gq_array_01[l,,] <- safe_ispline(input, derivs = 1)
  }
  
  # Construct hazard design matrix for GQ points, uncensored and in treatment group
  des_gq_array_11 <- array(NA_real_, dim = c(n_nodes, length(t11), df))
  for (l in seq_len(n_nodes)) {
    input <- t11/2 * (nodes[l] + 1)
    des_gq_array_11[l,,] <- safe_ispline(input, derivs = 1)
  }
  
  # Do the same for historical data (will create safe empty 3D arrays if length is 0)
  des_gq_array_01_hist <- array(0, dim = c(n_nodes, length(t01_hist), df))
  if (length(t01_hist) > 0) {
    for (l in seq_len(n_nodes)) {
      input <- t01_hist/2 * (nodes[l] + 1)
      des_gq_array_01_hist[l,,] <- safe_ispline(input, derivs = 1)
    }
  }
  
  des_gq_array_11_hist <- array(0, dim = c(n_nodes, length(t11_hist), df))
  if (length(t11_hist) > 0) {
    for (l in seq_len(n_nodes)) {
      input <- t11_hist/2 * (nodes[l] + 1)
      des_gq_array_11_hist[l,,] <- safe_ispline(input, derivs = 1)
    }
  }
  
  # Create objects for RMST
  t_nodes <- (t_upper / 2) * (nodes + 1)
  t_upper_design <- safe_ispline(t_nodes, derivs = 0)
  
  des_gq_array_rmst <- array(NA_real_, dim = c(n_nodes, n_nodes, df))
  for (l in seq_len(n_nodes)) {
    input <- t_nodes/2 * (nodes[l] + 1)
    des_gq_array_rmst[l,,] <- safe_ispline(input, derivs = 1)
  }
  
  data_list <- list(
    'n'         = n
    , 'n00'     = length(t00)
    , 'n01'     = length(t01)
    , 'n10'     = length(t10)
    , 'n11'     = length(t11)
    , 'df'      = df
    , 't00'     = t00
    , 't01'     = t01
    , 't10'     = t10
    , 't11'     = t11
    , 't_max'   = t_max
    , 't00_haz' = safe_ispline(t00, derivs = 1)
    , 't01_haz' = safe_ispline(t01, derivs = 1)
    , 't10_haz' = safe_ispline(t10, derivs = 1)
    , 't11_haz' = safe_ispline(t11, derivs = 1)
    , 't00_cum' = safe_ispline(t00, derivs = 0)
    , 't01_cum' = safe_ispline(t01, derivs = 0)
    , 't10_cum' = safe_ispline(t10, derivs = 0)
    , 't11_cum' = safe_ispline(t11, derivs = 0)
    , 'n_nodes' = n_nodes
    , 'nodes'   = nodes
    , 'weights' = weights
    , 'des_gq_01'  = des_gq_array_01
    , 'des_gq_11'  = des_gq_array_11
    , 'fix_tau' = ifelse(length(tau_val) == 1, 1, 0)
    , 'tau_val' = as.array(tau_val)
    , 'fix_epsilon' = ifelse(length(epsilon_val) == 1, 1, 0)
    , 'epsilon_val' = as.array(epsilon_val)
    , 'sd_epsilon' = sd_epsilon
    , 'rho'     = rho
    , 'sd_beta' = sd_beta
    , 'limit_tau' = ifelse(length(tau_max) == 1, 1, 0)
    , 'tau_max' = as.array(tau_max)
    , 'prior_tau' = prior_tau_int      # Passing the mapped integer to Stan
    , 'mixture_weight' = mixture_weight
    , 'sigma_tau' = sigma_tau
    , 't_upper' = t_upper
    , 't_nodes' = t_nodes
    , 't_upper_design' = t_upper_design
    , 'des_gq_rmst' = des_gq_array_rmst
    
    # Historical data (will be safe 0-length objects if a0=0 or omitted)
    , 'n_hist'  = n_hist
    , 'n00_hist'  = length(t00_hist)      
    , 'n01_hist'  = length(t01_hist)      
    , 'n10_hist'  = length(t10_hist)      
    , 'n11_hist'  = length(t11_hist)             
    , 't00_hist' = t00_hist         
    , 't01_hist' = t01_hist             
    , 't10_hist' = t10_hist                    
    , 't11_hist' = t11_hist                 
    , 't00_haz_hist' = safe_ispline(t00_hist, derivs = 1)
    , 't01_haz_hist' = safe_ispline(t01_hist, derivs = 1)
    , 't10_haz_hist' = safe_ispline(t10_hist, derivs = 1)
    , 't11_haz_hist' = safe_ispline(t11_hist, derivs = 1)
    , 't00_cum_hist' = safe_ispline(t00_hist, derivs = 0)
    , 't01_cum_hist' = safe_ispline(t01_hist, derivs = 0)
    , 't10_cum_hist' = safe_ispline(t10_hist, derivs = 0)
    , 't11_cum_hist' = safe_ispline(t11_hist, derivs = 0)
    , 'des_gq_01_hist' = des_gq_array_01_hist
    , 'des_gq_11_hist' = des_gq_array_11_hist  
    , 'a0' = a0
  )
  
  return(data_list)
}

#' Sample from posterior distribution of Bayesian cubic spline model
#' 
#' This is a wrapper function to obtain posterior samples for Bayesian cubic spline model
#' 
#' @param t vector of observed failure times
#' @param delta vector of event indicators
#' @param z vector of treatment indicators
#' @param t_hist vector of observed failure times in historical data (Optional)
#' @param delta_hist vector of event indicators in historical data (Optional)
#' @param z_hist vector of treatment indicators in historical data (Optional)
#' @param knots knot locations
#' @param a0 discounting parameter
#' @param prior_tau character: "uniform" (default), "mixture", or "half-normal" prior on tau
#' @param mixture_weight numeric: probability weight for the point mass spike at 0 (used if prior_tau == "mixture")
#' @param sigma_tau numeric: standard deviation for half-normal prior (used if prior_tau == "half-normal")
#' @param ... other arguments to pass onto `cmdstanr::sample`
#' 
#' @return an object of class `draws_df` giving posterior draws. See the `posterior` package for more details.
#' 
bayes_dte.mcmc <- function(
    t, delta, z, 
    t_hist = numeric(0), delta_hist = numeric(0), z_hist = numeric(0), 
    knots, n_nodes, 
    tau_val = numeric(0), tau_max = numeric(0), epsilon_val = numeric(0), 
    sd_epsilon = 10, rho = 0, sd_beta = 10, 
    t_upper, a0 = 0, 
    prior_tau = c("uniform", "mixture", "half-normal"), 
    mixture_weight = 0.5, sigma_tau = 1.0, ...) {
  
  # Get standata to pass to model (prior_tau is evaluated inside get.standata.dte)
  standata <- get.standata.dte(t, delta, z, t_hist, delta_hist, z_hist, 
                               knots, n_nodes, 
                               tau_val, tau_max, epsilon_val, 
                               sd_epsilon, rho, sd_beta, 
                               t_upper, a0,
                               prior_tau, mixture_weight, sigma_tau)
  
  # Obtain posterior draws
  fit_dte <- bayes.dte$sample(
    data = standata,
    # seed = 123,
    # chains = 1,
    # parallel_chains = 1,
    # iter_warmup = 2000,
    # iter_sampling = 10000,
    ...
  )
  
  return(fit_dte)
}