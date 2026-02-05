// Log hazard function for control group
// This function is defined by the derivative of an I-spline, aka an M-spline.
// The coefficients must be positive
vector log_hazard_ctrl(matrix design_haz, vector coef) {
  return log(design_haz * coef);
}

vector log_hazard_ratio(vector t, real tau, real epsilon, real beta_trt) {
  int T = size(t);
  vector[T] log_haz_ratio = rep_vector(0, T);

  // Precompute constants
  vector[4] alpha;
  alpha[1] = (pow(epsilon, 3) - 3 * epsilon * square(tau) + exp(beta_trt) * square(tau) * (3 * epsilon + 2 * tau) - 2 * tau^3);
  alpha[2] = (6 * tau * (epsilon - exp(beta_trt) * (epsilon + tau) + tau));
  alpha[3] = (3 * (-epsilon + exp(beta_trt) * (epsilon + 2 * tau) - 2 * tau));
  alpha[4] = (2 * (1 - exp(beta_trt)));

  for (j in 1:T) {
    if (t[j] > tau && t[j] <= tau + epsilon) {
      log_haz_ratio[j] += log((alpha[1] + alpha[2] * t[j] + alpha[3] * t[j]^2 + alpha[4] * t[j]^3) / pow(epsilon, 3));
    } else if (t[j] > tau + epsilon) {
      log_haz_ratio[j] += beta_trt;
    }
  }

  return log_haz_ratio;
}

// This is the quintic polynomial version of the hazard ratio
// vector log_hazard_ratio(vector t, real tau, real epsilon, real beta_trt) {
//   int T = size(t);
//   vector[T] log_haz_ratio = rep_vector(0, T);
//   
//   // Precompute constants
//   vector[6] alpha;
//   real delta = exp(beta_trt) - 1;
//   alpha[1] = 1 - delta * (10 * pow(tau, 3) / pow(epsilon, 3) + 15 * pow(tau, 4) / pow(epsilon, 4) + 6 * pow(tau, 5) / pow(epsilon, 5));
//   alpha[2] = 30 * square(tau) * delta * pow(tau + epsilon, 2) / pow(epsilon, 5);
//   alpha[3] = -30 * tau * delta * (tau + epsilon) * (2 * tau + epsilon) / pow(epsilon, 5);
//   alpha[4] = 10 * delta * (square(epsilon) + 6 * epsilon * tau + 6 * square(tau)) / pow(epsilon, 5);
//   alpha[5] = -15 * delta * (2 * tau + epsilon) / pow(epsilon, 5);
//   alpha[6] = 6 * delta / pow(epsilon, 5);
//   
//   for (j in 1:T) {
//     if (t[j] > tau && t[j] <= tau + epsilon) {
//       log_haz_ratio[j] += log((alpha[1] + alpha[2] * t[j] + alpha[3] * t[j]^2 + alpha[4] * t[j]^3 + alpha[5] * t[j]^4 + alpha[6] * t[j]^5));
//     } else if (t[j] > tau + epsilon) {
//       log_haz_ratio[j] += beta_trt;
//     }
//   }
// 
//   return log_haz_ratio;
// }


// Cumulative hazard function for the control group
// This function is defined by an I-spline
// The coefficients must be positive
vector cumulative_hazard_ctrl(matrix design_cum, vector coef) {
   return design_cum * coef;
}

// Return indices i such that x[i] <= thr
array[] int indices_le(vector x, real thr) {
  int N = size(x);
  array[N] int tmp;
  int n = 0;

  for (i in 1:N) {
    if (x[i] <= thr) {
      n += 1;
      tmp[n] = i;
    }
  }

  array[n == 0 ? 0 : n] int idx;
  if (n > 0) {
    for (i in 1:n) idx[i] = tmp[i];
  }
  return idx;
}

// Return indices i such that x[i] > thr
array[] int indices_gt(vector x, real thr) {
  int N = size(x);
  array[N] int tmp;
  int n = 0;

  for (i in 1:N) {
    if (x[i] > thr) {
      n += 1;
      tmp[n] = i;
    }
  }

  array[n == 0 ? 0 : n] int idx;
  if (n > 0) {
    for (i in 1:n) idx[i] = tmp[i];
  }
  return idx;
}


// Cumulative hazard function for the treatment group
// Integral computed using gaussian quadrature, with nodes and weights provided by user in R wrapper
// design cum is the full cumulative hazard design matrix corresponding to the vector t
vector cumulative_hazard_trt(matrix design_cum, vector t,
                             vector coef, real tau, real epsilon, real beta_trt,
                             array[,,] real des_gq,
                             // array[] matrix des_gq,
                             data vector nodes, data vector weights) {
  // Dimensions
  int n = size(t);
  int df = cols(design_cum);
  int L = size(nodes);
  
  // Initialize cumulative hazard vector
  vector[n] cum_haz_vec;
  
  // Get indices of observations less than and greater than tau
  int n_le = size(indices_le(t, tau));
  int n_gt = size(indices_gt(t, tau));
  array[n_le] int idx_le = indices_le(t, tau);  // indices with t <= tau
  array[n_gt] int idx_gt = indices_gt(t, tau);  // indices with t > tau
  
  // Set values for times before tau to baseline values
  cum_haz_vec[idx_le] = cumulative_hazard_ctrl(design_cum[idx_le,], coef);
  
  if (n_gt == 0) {
    return cum_haz_vec;
  }
  
  // Initialize hazard matrix
  matrix[n_gt, L] haz_mat;

  for (i in 1:L) {
    // Define submatrix used for GQ
    matrix[n_gt, df] des_sub;
    for (k in 1:n_gt) {
      // des_sub[k] = des_gq[i][idx_gt[k],];
      des_sub[k] = to_row_vector(des_gq[idx_gt[k], ,i]);
    }
    
    // log_hazard_ctrl(des_sub, coef) should return a vector of length n_gt
    vector[n_gt] log_haz_ctrl_vec = log_hazard_ctrl(des_sub, coef);
    vector[n_gt] trt_time = 0.5 * t[idx_gt] .* (nodes[i] + 1);
    for (k in 1:n_gt) {
      if (trt_time[k] > t[idx_gt][k])
        trt_time[k] = t[idx_gt][k];
    }
    vector[n_gt] log_haz_ratio_vec = log_hazard_ratio(trt_time, tau, epsilon, beta_trt);
    
    // combine and exponentiate
    haz_mat[,i] = exp(log_haz_ctrl_vec + log_haz_ratio_vec);
  }
  
  // Set values for times greater than tau
  cum_haz_vec[idx_gt] = (t[idx_gt] / 2) .* (haz_mat * weights);
  
  return cum_haz_vec;
}

// RMST
// This function computes the RMST in the control group
// The integral is computed from 0 to t_upper
real rmst_ctrl(matrix t_upper_design, real t_upper, vector coef,
               data vector nodes, data vector weights) {
  // t_upper_design: N x df matrix, row i is the basis evaluated at
  // u_i = t_upper/2 * (1 + nodes[i])
  int N = size(nodes);
  int df = cols(t_upper_design); // number of basis columns (should match length(coef))

  vector[N] surv_vec;

  for (i in 1:N) {
    row_vector[df] des_row = t_upper_design[i]; // row i => basis at node i
    real cum_haz = des_row * coef;              // scalar cumulative hazard
    surv_vec[i] = (t_upper / 2) * exp(-cum_haz); // quadrature integrand * Jacobian
  }

  return dot_product(surv_vec, weights); // weighted sum -> quadrature approx of integral
}

// This function computes the RMST in the treatment group
// The integral is computed from 0 to t_upper
// RMST for treatment group using cumulative_hazard_trt()
real rmst_trt(matrix t_upper_design_cum, vector t_nodes, real t_upper,
              vector coef, real tau, real epsilon, real beta_trt,
              array[,,] real des_gq,
              data vector nodes, data vector weights) {
                
  int N = size(nodes);

  // compute cumulative hazard at each quadrature node for the treatment group
  vector[N] cum_haz_trt = cumulative_hazard_trt(t_upper_design_cum, t_nodes,
                                                coef, tau, epsilon, beta_trt,
                                                des_gq, nodes, weights);

  // survival values at nodes
  vector[N] surv = exp(-cum_haz_trt);

  // quadrature weighted sum -> RMST
  return (t_upper / 2) * dot_product(surv, weights);
}
