// THIS IS THE MOST UP TO DATE VERSION OF THE STAN FILE THAT SHOULD BE USED
functions {
 #include dte_funs.stan
}
data {
  int<lower=0>                    n;              // total number of observations
  int<lower=0>                    n00;            // number of observations censored and in control group
  int<lower=0>                    n01;            // number of observations censored and in treatment group
  int<lower=0>                    n10;            // number of observations uncensored and in control group
  int<lower=0>                    n11;            // number of observations uncensored and in treatment group
  int<lower=4>                    df;             // degrees of freedom to use for cumulative hazard spline
  vector<lower=0>[n00]            t00;            // observed failure times censored and in control group
  vector<lower=0>[n01]            t01;            // observed failure times censored and in treatment group
  vector<lower=0>[n10]            t10;            // observed failure times uncensored and in control group
  vector<lower=0>[n11]            t11;            // observed failure times uncensored and in treatment group
  real<lower=0>                   t_max;          // maximum observed failure time to use for uniform distribution on tau
  matrix<lower=0>[n00, df]        t00_haz;        // design matrix of observed failure times censored and in control group
  matrix<lower=0>[n01, df]        t01_haz;        // design matrix of observed failure times censored and in treatment group
  matrix<lower=0>[n10, df]        t10_haz;        // design matrix of observed failure times uncensored and in control group
  matrix<lower=0>[n11, df]        t11_haz;        // design matrix of observed failure times uncensored and in treatment group
  matrix<lower=0>[n00, df]        t00_cum;        // design matrix of observed failure times censored and in control group
  matrix<lower=0>[n01, df]        t01_cum;        // design matrix of observed failure times censored and in treatment group
  matrix<lower=0>[n10, df]        t10_cum;        // design matrix of observed failure times uncensored and in control group
  matrix<lower=0>[n11, df]        t11_cum;        // design matrix of observed failure times uncensored and in treatment group
  int<lower=1>                    n_nodes;        // number of nodes to use for Gaussian Quadrature (GQ)
  vector[n_nodes]                 nodes;          // node locations for GQ
  vector[n_nodes]                 weights;        // weights for GQ
  array[n01, df, n_nodes] real<lower=0> des_gq_01;// array of baseline hazard design matrices for GQ (rows: subject, columns: df for M-spline, shelves: nodes of GQ), censored and in treatment group   
  array[n11, df, n_nodes] real<lower=0> des_gq_11;// array of baseline hazard design matrices for GQ (rows: subject, columns: df for M-spline, shelves: nodes of GQ), uncensored and in treatment group
  int<lower=0,upper=1>            fix_tau;        // indicator for fixing tau (e.g. to set it to 0)
  array[fix_tau ? 1 : 0] real<lower=0> tau_val;   // value for fixed tau (only provided if fix_tau is 1)
  int<lower=0,upper=1>            fix_epsilon;    // indicator for fixing epsilon
  array[fix_epsilon ? 1 : 0] real<lower=0> epsilon_val; // value for fixed epsilon (only provided if fix_epsilon is 1)
  // real<lower=0>                   sd_tau;         // standard deviation for tau's normal prior distribution
  real<lower=0>                   sd_epsilon;     // standard deviation for epsilons's normal prior distribution
  real<lower=0>                   rho;            // correlation for AR(1) structure in baseline cum_haz coefficients
  real<lower=0>                   sd_beta;        // prior sd for baseline cum_haz coefficients
  int<lower=0,upper=1>            limit_tau;      // indicator for whether to use alternative uniform prior for tau
  array[limit_tau ? 1 : 0] real<lower=0> tau_max; // upper value for uniform(0, tau_max) prior
  real<lower=0>                   t_upper;        // upper limit of integration for computing rmst
  vector[n_nodes]                 t_nodes;        // vector equal to t_upper/2 * (nodes + 1)
  matrix[n_nodes, df]             t_upper_design; // cumulative hazard design matrix of t_nodes
  array[n_nodes, df, n_nodes] real<lower=0> des_gq_rmst; // array of baseline hazard design matrices for GQ (rows: subject, columns: df for M-spline, shelves: nodes of GQ), censored and in treatment group   
}
parameters {
  vector[df] log_coef;                                       // spline coefficients on log scale
  // array[fix_tau ? 0 : 1] real<lower=0>       tau_raw;        // treatment effect delay, optionally fixed
  real<lower=0, upper=1> tau_unit; // scaled between 0 and 1
  array[fix_epsilon ? 0 : 1] real<lower=0>   epsilon_raw;    // transition window, optionally fixed
  real beta_trt;                                             // log of post-delay treatment effect
}
transformed parameters {
    real<lower=0> tau;
    if (fix_tau) {
        tau = tau_val[1];
    } else if (limit_tau) {
        tau = tau_unit * tau_max[1];
    } else {
        tau = tau_unit * t_max;
    }
    
    real<lower=0> epsilon;
    if (fix_epsilon) {
        epsilon = epsilon_val[1];       // Use fixed value
    } else {
        epsilon = epsilon_raw[1];       // Use estimated parameter
    }
    
    vector[df] coef = exp(log_coef);

     // AR(1) covariance matrix for baseline cum_haz coefficients
    cov_matrix[df] Sigma;
    for (i in 1:df) {
      for (j in 1:df) {
        Sigma[i,j] = square(sd_beta) * pow(rho, abs(i - j));
      }
    }
    
    real rmst_ctrl_val = rmst_ctrl(t_upper_design, t_upper, coef, nodes, weights);
    real rmst_trt_val = rmst_trt(t_upper_design, t_nodes, t_upper,
                                 coef, tau, epsilon, beta_trt,
                                 des_gq_rmst, nodes, weights);
    real rmst_diff = rmst_trt_val - rmst_ctrl_val;
}
model {
  // AR(1) multivariate normal on baseline hazard coefficients to induce smoothness
  vector[df] mu0 = rep_vector(0, df);
  log_coef ~ multi_normal(mu0, Sigma);
  
  // Noninformative prior on treatment effect
  beta_trt ~ normal(0, 10);
  
  // Set prior for tau depending on whether it's fixed
  if (!fix_tau) {
    tau_unit ~ uniform(0, 1);
  }
  
  // Set prior for epsilon depending on whether it's fixed
  if (!fix_epsilon) {
    epsilon_raw ~ normal(0, sd_epsilon);
  }
  
  // Likelihood components //
  // Censored and in control group
  target += -cumulative_hazard_ctrl(t00_cum, coef);
  // Censored and in treatment group
  target += -cumulative_hazard_trt(t01_cum, t01, coef, tau, epsilon, beta_trt, des_gq_01, nodes, weights);   
  // Uncensored and in control group
  target += log_hazard_ctrl(t10_haz, coef) -
    cumulative_hazard_ctrl(t10_cum, coef);
    // Uncensored and in treatment group
  target += log_hazard_ctrl(t11_haz, coef) +
    log_hazard_ratio(t11, tau, epsilon, beta_trt) - 
    cumulative_hazard_trt(t11_cum, t11, coef, tau, epsilon, beta_trt, des_gq_11, nodes, weights);
}
generated quantities {
  real hr_post = exp(beta_trt);
}
