# Load libraries
library(survival)
library(ggplot2)
library(ggsurvfit)
library(dplyr)
library(bayesplot)
library(survminer)
library(posterior)
library(survRM2)
library(nph)
library(patchwork)
library(ggtext)

# Source Stan model and R wrappers
source("Stan-R/dte_wrappers_pp.R")

# Source RMST R functions
source("Stan-R/rmst_funs.R")
source("Stan-R/RMST_Cox.R")

# Source files to generate data
source("Stan-R/generate_data.R")


# Set parameters to use to generate the data
# Sample size
n <- 250

# Post-delay hazard ratio and delay parameters
hr_post <- 0.7
tau <- 4
epsilon <- 2

# Knot locations
knots <- c(6.55, 9.63, 15.63)
bdry_knots <- c(0, 47.47)

# Baseline cumulative hazard coefficients
coef <- c(0.063397835, 0.007338422, 0.289976930, 1.469330398, 0.842779259, 0.015180263, 0.011023538)

# Compute true difference in RMST at 36 months
rmst_ctrl <- rmst_ctrl_r(36, knots, bdry_knots, coef) # 13.06
rmst_trt <- rmst_trt_r(36, knots, bdry_knots, coef, tau, epsilon, hr_post) # 15.73
rmst_diff <- rmst_trt - rmst_ctrl # 2.68

# Set seed to generate the data
set.seed(125)

# Ensure that max survival time in each group is greater than 36
t_max_min <- 0
while (t_max_min <= 36) {
  # Generate survival times
  t_control <- suppressWarnings(generate_data(n, knots, bdry_knots, coef))
  
  if (hr_post == 1) {
    t_trt <- suppressWarnings(generate_data(n, knots, bdry_knots, coef))
  } else{
    t_trt <- suppressWarnings(
      generate_data(n, knots, bdry_knots, coef, 
                    treatment = TRUE, tau, epsilon, hr_post)
    )
  }
  
  # Generate censoring times (e.g., from a uniform distribution)
  censoring_times <- runif(2*n, min = 0, max = 48)
  
  # Observed time is the minimum of event or censoring time
  observed_times <- pmin(c(t_control, t_trt), censoring_times)
  
  # Censoring indicator (1 = event occurred, 0 = censored)
  event <- ifelse(c(t_control, t_trt) <= censoring_times, 1, 0)
  
  # Create a survival data frame
  gendata <- data.frame(
    time = observed_times,
    event = event,
    group = c(rep("Control", n), rep("Treatment", n))
  ) 
  
  t_max_control <- max(gendata$time[which(gendata$group == "Control")])
  t_max_trt <- max(gendata$time[which(gendata$group == "Treatment")])
  t_max_min <- min(t_max_control, t_max_trt)
}


# Plot KM curves of generated data
km_fit <- survfit(Surv(time, event) ~ group, data = gendata)
ggsurvplot(km_fit,
           pval = TRUE,        # Add p-value for log-rank test
           risk.table = TRUE,  # Add risk table
           conf.int = TRUE,    # Add confidence intervals
           legend.title = "Group",
           xlab = "Time", ylab = "Survival Probability",
           title = "Kaplan-Meier Curves by Group")

# Set knot locations for fitting BayesDTE
knots.fit <- quantile(gendata %>% filter(group == "Control", event == 1) %>% pull(time),
                      probs = c(0.25, 0.5, 0.75))

# Provide initial values for sampler
initf <- function() list(log_coef = log(c(0.063397835, 0.007338422, 0.289976930, 1.469330398, 0.842779259, 0.015180263, 0.011023538)),
                         tau = 4,
                         # epsilon = 2,
                         beta_trt = 0)

# Fit BayesDTE with epsilon fixed at 2 and Uniform(0,12) prior on tau
fit.dte <- bayes_dte.mcmc(gendata$time, gendata$event, ifelse(gendata$group == "Treatment", 1, 0),
                          knots = knots.fit, n_nodes = 20, 
                          epsilon_val = 2, 
                          prior_tau = "uniform", tau_max = 12, 
                          rho = 0.5, sd_beta = 10, t_upper = 36,
                          chains = 4, parallel_chains = 4,
                          iter_warmup = 2000, iter_sampling = 2500,
                          refresh = 25, init = initf)


# Summary
summary.dte <- summary(as_draws_df(fit.dte), 
                       mean, sd, ~quantile(.x, probs = c(0.025, 0.975)), rhat, ess_bulk, ess_tail) 

summary.dte %>% 
  filter(variable %in% c("hr_post", "tau", "rmst_ctrl_val", "rmst_trt_val", "rmst_diff")) %>% 
  mutate(truth = c(tau, rmst_ctrl, rmst_trt, rmst_diff, hr_post))


# Posterior draws of the parameters that define S(t)
post <- as_tibble(as_draws_df(fit.dte$draws(c("coef", "tau", "epsilon", "beta_trt", "rmst_diff"))))
set.seed(123)
n_use <- nrow(post)                                      # nrow(post) to use every draw
post  <- slice_sample(post, n = n_use)
coef_draws <- as.matrix(post[, grep("^coef\\[", names(post))])   # S x df

# Spline basis used in the fit: knots.fit, with boundary knots c(0, t_max) where
# t_max = max(c(t, t_upper)), as computed in get.standata.dte()
fit_bdry_knots <- c(0, max(c(gendata$time, 36)))

# Time grid and Gaussian quadrature rule (same n_nodes as the fit)
t_grid  <- seq(0, fit_bdry_knots[2], length.out = 200)
n_nodes <- 20
gq      <- gauss.quad(n = n_nodes, kind = "legendre")

# Basis matrices for evaluating S at a vector of times: I-spline at the times
# (cum. haz.) and M-spline at the GQ points times/2 * (nodes + 1) (hazard)
surv_basis <- function(times) {
  u <- outer(times / 2, gq$nodes + 1)                                              # T x L
  list(times = times,
       I = iSpline(times, knots = knots.fit, Boundary.knots = fit_bdry_knots),     # T x df
       u = u,
       M = iSpline(as.vector(u), knots = knots.fit, Boundary.knots = fit_bdry_knots,
                   derivs = 1))                                                    # (T*L) x df
}
basis_grid <- surv_basis(t_grid)

# Hazard ratio, same cubic transition as hazard_ratio() in dte_funs.stan
hazard_ratio <- function(u, tau, epsilon, beta_trt) {
  eb <- exp(beta_trt)
  a1 <- epsilon^3 - 3 * epsilon * tau^2 + eb * tau^2 * (3 * epsilon + 2 * tau) - 2 * tau^3
  a2 <- 6 * tau * (epsilon - eb * (epsilon + tau) + tau)
  a3 <- 3 * (-epsilon + eb * (epsilon + 2 * tau) - 2 * tau)
  a4 <- 2 * (1 - eb)
  hr <- rep(1, length(u))
  mid <- u > tau & u <= tau + epsilon
  hr[mid] <- (a1 + a2 * u[mid] + a3 * u[mid]^2 + a4 * u[mid]^3) / epsilon^3
  hr[u > tau + epsilon] <- eb
  hr
}

# Treatment-arm survival at basis$times for one draw (mirrors cumulative_hazard_trt()):
# exact baseline cumulative hazard plus GQ of the treatment difference h0(u) * (HR(u) - 1),
# which is zero for u <= tau
surv_trt_draw <- function(coef, tau, epsilon, beta_trt, basis = basis_grid) {
  tt   <- basis$times
  H0   <- as.vector(basis$I %*% coef)                             # T
  h0   <- matrix(basis$M %*% coef, nrow = length(tt))             # T x L
  hr   <- matrix(hazard_ratio(as.vector(basis$u), tau, epsilon, beta_trt),
                 nrow = length(tt))
  H    <- H0 + (tt / 2) * as.vector((h0 * (hr - 1)) %*% gq$weights)
  exp(-H)
}

# T x S matrices of posterior survival curves for each arm
surv_draws <- function(basis) {
  list(ctrl = exp(-basis$I %*% t(coef_draws)),
       trt  = sapply(seq_len(n_use), function(s)
         surv_trt_draw(coef_draws[s, ], post$tau[s], post$epsilon[s], post$beta_trt[s], basis)))
}
S_grid <- surv_draws(basis_grid)
S_ctrl <- S_grid$ctrl
S_trt  <- S_grid$trt

# Pointwise and simultaneous (rank-envelope) credible bands
cred_bands <- function(S_mat, level = 0.95, time = t_grid) {
  n_s <- ncol(S_mat)
  pw  <- apply(S_mat, 1, quantile, probs = c((1 - level) / 2, (1 + level) / 2))
  # Rank each draw at each time; a draw's extremeness is its most extreme rank over t
  r       <- t(apply(S_mat, 1, rank, ties.method = "average"))     # G x S
  extreme <- apply(pmax(r, n_s + 1 - r), 2, max)
  k       <- ceiling(quantile(extreme, level, type = 1))
  S_sort  <- t(apply(S_mat, 1, sort))                                # G x S
  tibble(time      = time,
         mean      = rowMeans(S_mat),
         pw_lower  = pw[1, ], pw_upper  = pw[2, ],
         sim_lower = S_sort[, n_s + 1 - k], sim_upper = S_sort[, k],
         sim_coverage = mean(extreme <= k))
}

bands <- bind_rows(
  cred_bands(S_trt)  %>% mutate(trt_group = "Durvalumab + EP"),
  cred_bands(S_ctrl) %>% mutate(trt_group = "EP")
)

# Should be ~0.95 for each arm
distinct(bands, trt_group, sim_coverage)

# Posterior of S_trt(t) - S_ctrl(t) (columns of S_trt and S_ctrl are the same draws)
diff_bands <- cred_bands(S_trt - S_ctrl)

# RMST difference by truncation time L, by Gaussian quadrature as in rmst_trt()/rmst_ctrl():
# RMST(L) = L/2 * sum_l w_l S(L/2 * (node_l + 1)), with S at each node from the inner GQ rule
rmst_diff_draws <- function(L) {
  S_nodes <- surv_draws(surv_basis(as.vector(outer(L / 2, gq$nodes + 1))))       # (length(L)*n_nodes) x S
  A <- kronecker(t(gq$weights), diag(L / 2, nrow = length(L)))                    # length(L) x (length(L)*n_nodes)
  A %*% (S_nodes$trt - S_nodes$ctrl)
}
rmst_bands <- cred_bands(rmst_diff_draws(t_grid))

# Check against Stan's rmst_diff at t_upper = 36, draw by draw (should be ~0)
max(abs(rmst_diff_draws(36) - post$rmst_diff))

# Shared styling: arm colors/linetypes and band labels
arm_cols  <- c("Durvalumab + EP" = "#2E86C1", "EP" = "#CB4335")
arm_lty   <- c("Durvalumab + EP" = "solid",   "EP" = "dashed")
band_labs <- c("Pointwise 95% CrI", "Simultaneous 95% CrI")
band_alpha <- scale_alpha_manual(name = NULL, values = setNames(c(0.30, 0.15), band_labs),
                                 breaks = band_labs,
                                 guide = guide_legend(override.aes = list(fill = "grey30")))

## ---- Plot KM estimate with fitted curves and credible bands ----
fit_ipd <- survfit(Surv(time, event) ~ group, data = gendata)
mysurvplot <- ggsurvplot(
  fit_ipd,
  conf.int = FALSE,
  ggtheme = theme_classic(),
  data = gendata,
  legend.labs = c("EP", "Durvalumab + EP"),     # strata order: group=Control, group=Treatment
  palette = c("#CB4335", "#2E86C1"),
  legend.title="",
  risk.table = TRUE,
  risk.table.pos = "out",
  risk.table.col = "strata",
  tables.theme = theme_void() + theme(axis.text.y = element_text(size = 12)),
  break.time.by = 3,
  break.y.by = 0.2,
  xlim = c(0, 51),
  xlab = "Time (months)",
  ylab = "Probability of Overall Survival"
)

# Legend key for each arm: solid KM line in the arm color (passed in via the key's fill)
# with the black fitted line (solid or dashed) drawn on top
draw_key_km_fit <- function(data, params, size) {
  km  <- draw_key_path(transform(data, colour = fill, linetype = "solid", linewidth = 2.5 * linewidth),
                       params, size)
  fit <- draw_key_path(data, params, size)
  grid::grobTree(km, fit)
}

# Fitted curves: posterior mean of S(t), solid for Durvalumab + EP and dashed for EP.
# The arm legend shows KM color plus fitted linetype; ggsurvplot's own color legend is dropped.
mysurvplot$plot <- mysurvplot$plot +
  geom_line(data = bands, aes(time, mean, linetype = trt_group),
            color = "black", linewidth = 0.85, inherit.aes = FALSE, key_glyph = draw_key_km_fit) +
  scale_linetype_manual(name = NULL, values = arm_lty,
                        guide = guide_legend(order = 1, override.aes = list(fill = arm_cols))) +
  guides(color = "none") +
  labs(x = "Time (months)", y = "Probability of Overall Survival") +
  theme(legend.key.width = unit(1.6, "cm"))

# Add bands, filled by arm color; the band legend uses grey keys at matching transparency
band_layers <- unlist(lapply(names(arm_cols), function(arm) {
  b <- filter(bands, trt_group == arm)
  list(geom_ribbon(data = b, aes(x = time, ymin = sim_lower, ymax = sim_upper, alpha = band_labs[2]),
                   fill = arm_cols[arm], inherit.aes = FALSE),
       geom_ribbon(data = b, aes(x = time, ymin = pw_lower, ymax = pw_upper, alpha = band_labs[1]),
                   fill = arm_cols[arm], inherit.aes = FALSE))
}), recursive = FALSE)

mysurvplot$plot <- mysurvplot$plot + band_layers + band_alpha
# Move the ribbons underneath the KM steps and fitted lines
n_new <- length(band_layers)
mysurvplot$plot$layers <- c(tail(mysurvplot$plot$layers, n_new),
                            head(mysurvplot$plot$layers, -n_new))

scale <- 1.8
text_size <- 9*scale
mysurvplot$plot <- mysurvplot$plot +
  theme(axis.text = element_text(size = 14, color = "black"),
        legend.text = element_text(size = text_size*0.8, color = "black"),
        axis.title = element_text(size = text_size),
        axis.title.y = element_text(margin = margin(r = 10)),   # gap between y-axis title and tick labels
        axis.line = element_line(color = "black"))

# Drop the stray strata legend next to the risk table
mysurvplot$table <- mysurvplot$table + theme(legend.position = "none")

# Start the x axis exactly at 0 (no left padding) so the y axis meets it at t = 0;
# same scale on the risk table so the columns stay aligned, and clip = "off" so the
# t = 0 risk counts and axis labels aren't cut off at the panel edge
x_scale <- scale_x_continuous(breaks = seq(0, 51, 3), expand = expansion(mult = c(0, 0.02)))
mysurvplot$plot  <- mysurvplot$plot  + x_scale + coord_cartesian(xlim = c(0, 51), clip = "off")
mysurvplot$table <- mysurvplot$table + x_scale + coord_cartesian(xlim = c(0, 51), clip = "off") +
  theme(axis.text.y = element_markdown(margin = margin(r = 20)))   # room for the t = 0 counts

pdf("Example/fitted_plot.pdf", height = scale*4, width = scale*7)
free(mysurvplot$plot, type = "label", side = "l") / mysurvplot$table + plot_layout(heights = c(8,2))
dev.off()
