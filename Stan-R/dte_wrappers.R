library(cmdstanr)
library(rstan)
library(posterior)
library(statmod)
library(splines2)

bayes.dte <- cmdstanr::cmdstan_model("Stan-R/dte.stan")

#' Get stan data for BayesDTE
#' 
#' Obtains stan data for BayesDTE
#' 
#' @param t vector of observed failure times
#' @param delta vector of event indicators
#' @param z vector of treatment indicators
#' @param knots knot locations
#' @param n_nodes number of nodes to use for Gaussian quadrature
#' @param tau_val
#' @param epsilon_val
#'
#' @return a list giving data to be passed onto `rstan` or `cmdstanr` models.
get.standata.dte <- function(t, delta, z, knots, n_nodes, tau_val = numeric(0), tau_max = numeric(0), epsilon_val = numeric(0), sd_epsilon = 10, rho = 0, sd_beta = 10, t_upper) {
  n <- length(t)
  df <- length(knots) + 4
  
  # Compute Gaussian quadrature
  gq <- gauss.quad(n = n_nodes, kind = "legendre")
  nodes <- gq$nodes
  weights <- gq$weights
  
  t00 <- t[delta == 0 & z == 0]
  t01 <- t[delta == 0 & z == 1]
  t10 <- t[delta == 1 & z == 0]
  t11 <- t[delta == 1 & z == 1]
  
  t_max <- max(t)
  bdry_knots <- c(0, t_max)
  
  # Construct hazard design matrix for GQ points, censored and in treatment group
  des_gq_01 <- list()
  for (l in 1:length(nodes)) {
    input <- t01/2 * (nodes[l] + 1)
    D_l <- iSpline(input, knots=knots, Boundary.knots = bdry_knots, derivs = 1)
    des_gq_01[[l]] <- as.matrix(D_l)
  }
  
  # Fix ordering and provide as array for stan, censored and in treatment group
  des_gq_array_01 <- array(NA_real_, dim = c(length(t01), df, n_nodes))
  for (l in seq_len(n_nodes)) {
    des_gq_array_01[,,l] <- des_gq_01[[l]]
  }
  
  # Construct hazard design matrix for GQ points, uncensored and in treatment group
  des_gq_11 <- list()
  for (l in 1:length(nodes)) {
    input <- t11/2 * (nodes[l] + 1)
    D_l <- iSpline(input, knots=knots, Boundary.knots = bdry_knots, derivs = 1)
    des_gq_11[[l]] <- as.matrix(D_l)
  }
  
  # Fix ordering and provide as array for stan, uncensored and in treatment group
  des_gq_array_11 <- array(NA_real_, dim = c(length(t11), df, n_nodes))
  for (l in seq_len(n_nodes)) {
    des_gq_array_11[,,l] <- des_gq_11[[l]]
  }
  
  # Create objects for RMST
  t_nodes <- (t_upper / 2) * (nodes + 1)
  t_upper_design <- iSpline(t_nodes, knots = knots, Boundary.knots = bdry_knots)
  
  des_gq_rmst <- list()
  for (l in 1:length(nodes)) {
    input <- t_nodes/2 * (nodes[l] + 1)
    D_l <- iSpline(input, knots=knots, Boundary.knots = bdry_knots, derivs = 1)
    des_gq_rmst[[l]] <- as.matrix(D_l)
  }
  
  des_gq_array_rmst <- array(NA_real_, dim = c(n_nodes, df, n_nodes))
  for (l in seq_len(n_nodes)) {
    des_gq_array_rmst[,,l] <- des_gq_rmst[[l]]
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
    , 't00_haz' = iSpline(t00, knots = knots, Boundary.knots = bdry_knots, derivs = 1)
    , 't01_haz' = iSpline(t01, knots = knots, Boundary.knots = bdry_knots, derivs = 1)
    , 't10_haz' = iSpline(t10, knots = knots, Boundary.knots = bdry_knots, derivs = 1)
    , 't11_haz' = iSpline(t11, knots = knots, Boundary.knots = bdry_knots, derivs = 1)
    , 't00_cum' = iSpline(t00, knots = knots, Boundary.knots = bdry_knots)
    , 't01_cum' = iSpline(t01, knots = knots, Boundary.knots = bdry_knots)
    , 't10_cum' = iSpline(t10, knots = knots, Boundary.knots = bdry_knots)
    , 't11_cum' = iSpline(t11, knots = knots, Boundary.knots = bdry_knots)
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
    , 't_upper' = t_upper
    , 't_nodes' = t_nodes
    , 't_upper_design' = t_upper_design
    , 'des_gq_rmst' = des_gq_array_rmst
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
#' @param knots knot locations
#' @param c 75th percentile of historical data to be used for log basis functions
#' @param ... other arguments to pass onto `cmdstanr::sample`
#' 
#' @return an object of class `draws_df` giving posterior draws. See the `posterior` package for more details.
#' 
bayes_dte.mcmc <- function(t, delta, z, knots, n_nodes, tau_val = numeric(0), tau_max = numeric(0), epsilon_val = numeric(0), sd_epsilon = 10, rho = 0, sd_beta = 10, t_upper, ...) {
  # Get standata to pass to model
  standata <- get.standata.dte(t, delta, z, knots, n_nodes, tau_val, tau_max, epsilon_val, sd_epsilon, rho, sd_beta, t_upper)
  
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
