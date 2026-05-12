// Log hazard function for control group
// This function is defined by the derivative of an I-spline, aka an M-spline.
// The coefficients must be positive
vector log_hazard_ctrl(matrix design_haz, vector coef) {
  return log(design_haz * coef);
}

// Hazard ratio function (Returns raw ratio, NOT log!)
// Evaluates using Horner's method to avoid expensive pow() and exp() calls in the loop
vector hazard_ratio(vector t, real tau, real epsilon, real beta_trt) {
  int T = size(t);
  vector[T] hr = rep_vector(1.0, T); // Default to 1 for t <= tau

  real exp_beta = exp(beta_trt);
  real eps3 = pow(epsilon, 3);
  real tau2 = square(tau);

  // Pre-divide by eps3 here to save division inside the loop
  real a1 = (eps3 - 3 * epsilon * tau2 + exp_beta * tau2 * (3 * epsilon + 2 * tau) - 2 * tau2 * tau) / eps3;
  real a2 = (6 * tau * (epsilon - exp_beta * (epsilon + tau) + tau)) / eps3;
  real a3 = (3 * (-epsilon + exp_beta * (epsilon + 2 * tau) - 2 * tau)) / eps3;
  real a4 = (2 * (1 - exp_beta)) / eps3;

  for (j in 1:T) {
    if (t[j] > tau + epsilon) {
      hr[j] = exp_beta;
    } else if (t[j] > tau) {
      // Horner's method for O(N) cubic evaluation
      hr[j] = ((a4 * t[j] + a3) * t[j] + a2) * t[j] + a1;
    }
  }

  return hr;
}

// Cumulative hazard function for the control group
// This function is defined by an I-spline
// The coefficients must be positive
vector cumulative_hazard_ctrl(matrix design_cum, vector coef) {
   return design_cum * coef;
}

// Helper to precompute GQ times (moved from inner loops to initialization)
matrix make_gq_times(vector t, vector nodes) {
  int N = size(t);
  int L = size(nodes);
  matrix[N, L] gq_times;
  
  for (i in 1:L) {
    gq_times[, i] = 0.5 * t * (nodes[i] + 1);
    for (j in 1:N) {
      if (gq_times[j, i] > t[j]) {
        gq_times[j, i] = t[j];
      }
    }
  }
  return gq_times;
}

// Cumulative hazard function for the treatment group
// Integral computed using gaussian quadrature
// Uses the zero-difference trick to avoid bounds checking and array resizing
vector cumulative_hazard_trt(matrix design_cum, vector t,
                             vector coef, real tau, real epsilon, real beta_trt,
                             array[] matrix des_gq,
                             matrix gq_times,
                             vector weights) {
  int n = size(t);
  int L = size(weights);
  
  // 1. Exact baseline cumulative hazard for everyone
  vector[n] cum_haz = cumulative_hazard_ctrl(design_cum, coef);
  
  // 2. Add the GQ approximation of ONLY the treatment difference
  matrix[n, L] diff_mat;

  for (i in 1:L) {
    // Single matrix-vector multiply
    vector[n] haz_0 = des_gq[i] * coef; 
    
    // Get HR for all nodes
    vector[n] hr = hazard_ratio(gq_times[, i], tau, epsilon, beta_trt);
    
    // Where t <= tau, hr=1, so (hr - 1.0) = 0, perfectly zeroing out the effect
    diff_mat[, i] = haz_0 .* (hr - 1.0);
  }
  
  // Integrate the difference and add to baseline
  cum_haz += (t / 2.0) .* (diff_mat * weights);
  
  return cum_haz;
}

// RMST for control group
real rmst_ctrl(matrix t_upper_design, real t_upper, vector coef,
               vector nodes, vector weights) {
  int N = size(nodes);
  int df = cols(t_upper_design); 

  vector[N] surv_vec;

  for (i in 1:N) {
    row_vector[df] des_row = t_upper_design[i]; 
    real cum_haz = des_row * coef;              
    surv_vec[i] = (t_upper / 2) * exp(-cum_haz); 
  }

  return dot_product(surv_vec, weights); 
}

// RMST for treatment group
real rmst_trt(matrix t_upper_design_cum, vector t_nodes, real t_upper,
              vector coef, real tau, real epsilon, real beta_trt,
              array[] matrix des_gq, matrix gq_times_rmst,
              vector nodes, vector weights) {
                
  int N = size(nodes);

  vector[N] cum_haz_trt = cumulative_hazard_trt(t_upper_design_cum, t_nodes,
                                                coef, tau, epsilon, beta_trt,
                                                des_gq, gq_times_rmst, weights);

  vector[N] surv = exp(-cum_haz_trt);

  return (t_upper / 2) * dot_product(surv, weights);
}