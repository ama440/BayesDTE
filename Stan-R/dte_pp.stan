// This is the most up to date version of the stan file
// It now includes an option to incorporate historical data using the power prior
// Optimized for HMC sampling via vectorization and zero-difference integration

functions {
  #include dte_funs.stan
}
data {
  int<lower=0>                    n;               // total number of observations
  int<lower=0>                    n00;             // number of observations censored and in control group
  int<lower=0>                    n01;             // number of observations censored and in treatment group
  int<lower=0>                    n10;             // number of observations uncensored and in control group
  int<lower=0>                    n11;             // number of observations uncensored and in treatment group
  int<lower=4>                    df;              // degrees of freedom to use for cumulative hazard spline
  vector<lower=0>[n00]            t00;             // observed failure times censored and in control group
  vector<lower=0>[n01]            t01;             // observed failure times censored and in treatment group
  vector<lower=0>[n10]            t10;             // observed failure times uncensored and in control group
  vector<lower=0>[n11]            t11;             // observed failure times uncensored and in treatment group
  real<lower=0>                   t_max;           // maximum observed failure time to use for uniform distribution on tau
  matrix<lower=0>[n00, df]        t00_haz;         // design matrix of observed failure times censored and in control group
  matrix<lower=0>[n01, df]        t01_haz;         // design matrix of observed failure times censored and in treatment group
  matrix<lower=0>[n10, df]        t10_haz;         // design matrix of observed failure times uncensored and in control group
  matrix<lower=0>[n11, df]        t11_haz;         // design matrix of observed failure times uncensored and in treatment group
  matrix<lower=0>[n00, df]        t00_cum;         // design matrix of observed failure times censored and in control group
  matrix<lower=0>[n01, df]        t01_cum;         // design matrix of observed failure times censored and in treatment group
  matrix<lower=0>[n10, df]        t10_cum;         // design matrix of observed failure times uncensored and in control group
  matrix<lower=0>[n11, df]        t11_cum;         // design matrix of observed failure times uncensored and in treatment group
  int<lower=1>                    n_nodes;         // number of nodes to use for Gaussian Quadrature (GQ)
  vector[n_nodes]                 nodes;           // node locations for GQ
  vector[n_nodes]                 weights;         // weights for GQ
  
  // NOTE: Transposed array shapes for speed!
  array[n_nodes] matrix<lower=0>[n01, df] des_gq_01; // array of baseline hazard design matrices for GQ (shelves: nodes, rows: subject, cols: df)
  array[n_nodes] matrix<lower=0>[n11, df] des_gq_11; // array of baseline hazard design matrices for GQ (shelves: nodes, rows: subject, cols: df)
  
  int<lower=0,upper=1>            fix_tau;         // indicator for fixing tau (e.g. to set it to 0)
  array[fix_tau ? 1 : 0] real<lower=0> tau_val;    // value for fixed tau (only provided if fix_tau is 1)
  int<lower=0,upper=1>            fix_epsilon;     // indicator for fixing epsilon
  array[fix_epsilon ? 1 : 0] real<lower=0> epsilon_val; // value for fixed epsilon (only provided if fix_epsilon is 1)
  real<lower=0>                   sd_epsilon;      // standard deviation for epsilons's normal prior distribution
  real<lower=0>                   rho;             // correlation for AR(1) structure in baseline cum_haz coefficients
  real<lower=0>                   sd_beta;         // prior sd for baseline cum_haz coefficients
  
  // --- tau ---
  int<lower=0,upper=1>            limit_tau;       // indicator for whether to use alternative uniform prior for tau
  array[limit_tau ? 1 : 0] real<lower=0> tau_max;  // upper value for uniform(0, tau_max) prior
  int<lower=0, upper=2>           prior_tau;       // 0: uniform, 1: mixture (spike-and-slab), 2: half-normal
  real<lower=0, upper=1>          mixture_weight;  // probability weight for the point mass spike at 0 (used if prior_tau == 1)
  real<lower=0>                   sigma_tau;       // standard deviation for half-normal (used if prior_tau == 2)
  
  // RMST
  real<lower=0>                   t_upper;         // upper limit of integration for computing rmst
  vector[n_nodes]                 t_nodes;         // vector equal to t_upper/2 * (nodes + 1)
  matrix[n_nodes, df]             t_upper_design;  // cumulative hazard design matrix of t_nodes
  
  // NOTE: Transposed array shapes for speed
  array[n_nodes] matrix<lower=0>[n_nodes, df] des_gq_rmst; 
  
  // Include historical data
  int<lower=0>                    n_hist;              
  int<lower=0>                    n00_hist;            
  int<lower=0>                    n01_hist;            
  int<lower=0>                    n10_hist;            
  int<lower=0>                    n11_hist;            
  vector<lower=0>[n00_hist]       t00_hist;            
  vector<lower=0>[n01_hist]       t01_hist;            
  vector<lower=0>[n10_hist]       t10_hist;            
  vector<lower=0>[n11_hist]       t11_hist;            
  matrix<lower=0>[n00_hist, df]   t00_haz_hist;        
  matrix<lower=0>[n01_hist, df]   t01_haz_hist;        
  matrix<lower=0>[n10_hist, df]   t10_haz_hist;        
  matrix<lower=0>[n11_hist, df]   t11_haz_hist;        
  matrix<lower=0>[n00_hist, df]   t00_cum_hist;        
  matrix<lower=0>[n01_hist, df]   t01_cum_hist;        
  matrix<lower=0>[n10_hist, df]   t10_cum_hist;        
  matrix<lower=0>[n11_hist, df]   t11_cum_hist;        

  // NOTE: Transposed array shapes for speed
  array[n_nodes] matrix<lower=0>[n01_hist, df] des_gq_01_hist;    
  array[n_nodes] matrix<lower=0>[n11_hist, df] des_gq_11_hist;    

  real<lower=0,upper=1> a0;          // discounting parameter for power prior
} 

transformed data {
  // Precompute AR(1) Cholesky covariance to save leapfrog iterations
  matrix[df, df] Sigma;
  for (i in 1:df) {
    for (j in 1:df) {
      Sigma[i,j] = square(sd_beta) * pow(rho, abs(i - j));
    }
  }
  matrix[df, df] L_Sigma = cholesky_decompose(Sigma);
  vector[df] mu0 = rep_vector(0, df);
  
  // Precompute numerical integration time points for all groups
  matrix[n01, n_nodes] gq_times_01 = make_gq_times(t01, nodes);
  matrix[n11, n_nodes] gq_times_11 = make_gq_times(t11, nodes);
  matrix[n_nodes, n_nodes] gq_times_rmst = make_gq_times(t_nodes, nodes);
  
  matrix[n01_hist, n_nodes] gq_times_01_hist = make_gq_times(t01_hist, nodes);
  matrix[n11_hist, n_nodes] gq_times_11_hist = make_gq_times(t11_hist, nodes);
}

parameters {
  vector[df] log_coef;                                       // spline coefficients on log scale
  real<lower=0, upper=1> tau_unit;                           // Used when prior_tau is 0 or 1 (scaled between 0 and 1)
  real<lower=0>          tau_hn;                             // Used when prior_tau is 2 (unbounded upper limit for true half-normal)
  array[fix_epsilon ? 0 : 1] real<lower=0>  epsilon_raw;     // transition window, optionally fixed
  real beta_trt;                                             // log of post-delay treatment effect
}

transformed parameters {
    real<lower=0> tau;
    if (fix_tau) {
        tau = tau_val[1];
    } else if (prior_tau == 2) {
        tau = tau_hn;
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
    
    real rmst_ctrl_val = rmst_ctrl(t_upper_design, t_upper, coef, nodes, weights);
    real rmst_trt_val = rmst_trt(t_upper_design, t_nodes, t_upper,
                                 coef, tau, epsilon, beta_trt,
                                 des_gq_rmst, gq_times_rmst, nodes, weights);
    real rmst_diff = rmst_trt_val - rmst_ctrl_val;
}

model {
  // Faster Cholesky AR(1) multivariate normal
  log_coef ~ multi_normal_cholesky(mu0, L_Sigma);
  
  // Noninformative prior on treatment effect
  beta_trt ~ normal(0, 10);
  
  // Set prior for tau depending on chosen option
  if (!fix_tau) {
    if (prior_tau == 0 || prior_tau == 1) {
      // For Option 0 (Uniform) and Option 1 (Marginalized Mixture), 
      // the parameter tau_unit just needs a standard uniform prior.
      tau_unit ~ uniform(0, 1);
      tau_hn ~ normal(0, 1); // Dummy prior
    } else if (prior_tau == 2) {
      // Option 2: Half-normal prior
      tau_hn ~ normal(0, sigma_tau);
      tau_unit ~ uniform(0, 1); // Dummy prior
    }
  }
  
  // Set prior for epsilon depending on whether it's fixed
  if (!fix_epsilon) {
    epsilon_raw ~ normal(0, sd_epsilon);
  }
  
  // Likelihood components for current data //
  
  // --- 1. Components INDEPENDENT of tau ---
  // We compute these once to save computation time
  
  target += -cumulative_hazard_ctrl(t00_cum, coef);
  target += log_hazard_ctrl(t10_haz, coef) - cumulative_hazard_ctrl(t10_cum, coef);
  target += log_hazard_ctrl(t11_haz, coef); // Baseline hazard for treatment uncensored
  
  if (a0 > 0) {
    target += -a0*cumulative_hazard_ctrl(t00_cum_hist, coef);
    target += a0*(log_hazard_ctrl(t10_haz_hist, coef) - cumulative_hazard_ctrl(t10_cum_hist, coef));
    target += a0*log_hazard_ctrl(t11_haz_hist, coef);
  }

  // --- 2. Components DEPENDENT on tau ---
  
  if (!fix_tau && prior_tau == 1) {
    // MARGINALIZED MIXTURE LOGIC
    real log_lik_null = 0;
    real log_lik_dte = 0;

    // A. Compute Null state (tau = 0)
    log_lik_null += sum(-cumulative_hazard_trt(t01_cum, t01, coef, 0.0, epsilon, beta_trt, des_gq_01, gq_times_01, weights));
    log_lik_null += sum(log(hazard_ratio(t11, 0.0, epsilon, beta_trt)) - cumulative_hazard_trt(t11_cum, t11, coef, 0.0, epsilon, beta_trt, des_gq_11, gq_times_11, weights));

    // B. Compute DTE state (tau = sampled parameter)
    log_lik_dte += sum(-cumulative_hazard_trt(t01_cum, t01, coef, tau, epsilon, beta_trt, des_gq_01, gq_times_01, weights));
    log_lik_dte += sum(log(hazard_ratio(t11, tau, epsilon, beta_trt)) - cumulative_hazard_trt(t11_cum, t11, coef, tau, epsilon, beta_trt, des_gq_11, gq_times_11, weights));

    // C. Add Historical Data (if applicable)
    if (a0 > 0) {
      log_lik_null += sum(-a0*cumulative_hazard_trt(t01_cum_hist, t01_hist, coef, 0.0, epsilon, beta_trt, des_gq_01_hist, gq_times_01_hist, weights));
      log_lik_null += sum(a0*(log(hazard_ratio(t11_hist, 0.0, epsilon, beta_trt)) - cumulative_hazard_trt(t11_cum_hist, t11_hist, coef, 0.0, epsilon, beta_trt, des_gq_11_hist, gq_times_11_hist, weights)));

      log_lik_dte += sum(-a0*cumulative_hazard_trt(t01_cum_hist, t01_hist, coef, tau, epsilon, beta_trt, des_gq_01_hist, gq_times_01_hist, weights));
      log_lik_dte += sum(a0*(log(hazard_ratio(t11_hist, tau, epsilon, beta_trt)) - cumulative_hazard_trt(t11_cum_hist, t11_hist, coef, tau, epsilon, beta_trt, des_gq_11_hist, gq_times_11_hist, weights)));
    }

    // D. Mix the likelihoods on the log scale
    target += log_mix(mixture_weight, log_lik_null, log_lik_dte);

  } else {
    // STANDARD LOGIC (For Uniform, Half-Normal, or Fixed tau)
    target += -cumulative_hazard_trt(t01_cum, t01, coef, tau, epsilon, beta_trt, des_gq_01, gq_times_01, weights);
    target += log(hazard_ratio(t11, tau, epsilon, beta_trt)) - cumulative_hazard_trt(t11_cum, t11, coef, tau, epsilon, beta_trt, des_gq_11, gq_times_11, weights);

    if (a0 > 0) {
      target += -a0*cumulative_hazard_trt(t01_cum_hist, t01_hist, coef, tau, epsilon, beta_trt, des_gq_01_hist, gq_times_01_hist, weights);
      target += a0*(log(hazard_ratio(t11_hist, tau, epsilon, beta_trt)) - cumulative_hazard_trt(t11_cum_hist, t11_hist, coef, tau, epsilon, beta_trt, des_gq_11_hist, gq_times_11_hist, weights));
    }
  }
}

generated quantities {
  real hr_post = exp(beta_trt);
  
  // --- NEW: Reconstructed Marginalized Quantities ---
  real prob_null;          // Posterior probability that tau = 0 for this draw
  int z_sim;               // Simulated discrete state (0 = Null/tau=0, 1 = DTE/tau>0)
  real rmst_diff_mix;      // The true RMST difference accounting for the mixture state
  
  if (!fix_tau && prior_tau == 1) {
    real log_lik_null_gq = 0;
    real log_lik_dte_gq = 0;
    
    // 1. Recompute the likelihoods for the saved parameter draw
    log_lik_null_gq += sum(-cumulative_hazard_trt(t01_cum, t01, coef, 0.0, epsilon, beta_trt, des_gq_01, gq_times_01, weights));
    log_lik_null_gq += sum(log(hazard_ratio(t11, 0.0, epsilon, beta_trt)) - cumulative_hazard_trt(t11_cum, t11, coef, 0.0, epsilon, beta_trt, des_gq_11, gq_times_11, weights));

    log_lik_dte_gq += sum(-cumulative_hazard_trt(t01_cum, t01, coef, tau, epsilon, beta_trt, des_gq_01, gq_times_01, weights));
    log_lik_dte_gq += sum(log(hazard_ratio(t11, tau, epsilon, beta_trt)) - cumulative_hazard_trt(t11_cum, t11, coef, tau, epsilon, beta_trt, des_gq_11, gq_times_11, weights));

    // Add historical data likelihoods if applicable
    if (a0 > 0) {
      log_lik_null_gq += sum(-a0*cumulative_hazard_trt(t01_cum_hist, t01_hist, coef, 0.0, epsilon, beta_trt, des_gq_01_hist, gq_times_01_hist, weights));
      log_lik_null_gq += sum(a0*(log(hazard_ratio(t11_hist, 0.0, epsilon, beta_trt)) - cumulative_hazard_trt(t11_cum_hist, t11_hist, coef, 0.0, epsilon, beta_trt, des_gq_11_hist, gq_times_11_hist, weights)));

      log_lik_dte_gq += sum(-a0*cumulative_hazard_trt(t01_cum_hist, t01_hist, coef, tau, epsilon, beta_trt, des_gq_01_hist, gq_times_01_hist, weights));
      log_lik_dte_gq += sum(a0*(log(hazard_ratio(t11_hist, tau, epsilon, beta_trt)) - cumulative_hazard_trt(t11_cum_hist, t11_hist, coef, tau, epsilon, beta_trt, des_gq_11_hist, gq_times_11_hist, weights)));
    }
    
    // 2. Compute posterior probabilities using Bayes' rule in log space
    // Using log space prevents floating point underflow errors
    real log_post_null = log(mixture_weight) + log_lik_null_gq;
    real log_post_dte = log(1.0 - mixture_weight) + log_lik_dte_gq;
    
    // Normalize to get the actual probability of the null state
    prob_null = exp(log_post_null - log_sum_exp(log_post_null, log_post_dte));
    
    // 3. Simulate the discrete state (bernoulli_rng takes the probability of success, i.e., state 1 / DTE)
    z_sim = bernoulli_rng(1.0 - prob_null);
    
    // 4. Calculate the final true RMST difference
    if (z_sim == 0) {
      // If the state is simulated as Null, calculate RMST using tau = 0
      real rmst_trt_val_null = rmst_trt(t_upper_design, t_nodes, t_upper,
                                        coef, 0.0, epsilon, beta_trt,
                                        des_gq_rmst, gq_times_rmst, nodes, weights);
      rmst_diff_mix = rmst_trt_val_null - rmst_ctrl_val; 
    } else {
      // If the state is simulated as DTE, use the continuous RMST calculated in transformed parameters
      rmst_diff_mix = rmst_diff; 
    }
    
  } else {
    // Standard pass-through for non-mixture models
    prob_null = 0.0; 
    z_sim = 1;
    rmst_diff_mix = rmst_diff;
  }
}
