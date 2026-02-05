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

# Source Stan model and R wrappers
source("Stan-R/dte_wrappers.R")

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
set.seed(123)

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

# Fit BayesDTE
fit.dte <- bayes_dte.mcmc(gendata$time, gendata$event, ifelse(gendata$group == "Treatment", 1, 0),
                          knots.fit, n_nodes = 10, epsilon_val = 2, tau_max = 12, rho = 0.5, sd_beta = 10, t_upper = 36,
                          chains = 4, parallel_chains = 4,
                          iter_warmup = 2000, iter_sampling = 2500,
                          refresh = 25, init = initf)


# Summary
summary.dte <- summary(as_draws_df(fit.dte), 
                       mean, sd, ~quantile(.x, probs = c(0.025, 0.975)), rhat, ess_bulk, ess_tail) 

summary.dte %>% 
  filter(variable %in% c("hr_post", "tau", "rmst_ctrl_val", "rmst_trt_val", "rmst_diff"))


# Visualize fitted curves; requires Rstan
library(rstan)
dte.stan  <- stan_model("Stan-R/RStan/test.stan")
functions <- expose_stan_functions(dte.stan)

# Get design matrix and posterior means of baseline hazard coefficients
design <- iSpline(gendata$t[which(gendata$group=="Control")], knots = knots.fit, Boundary.knots = c(0, max(gendata$t)))
idx <- grep("^coef", summary.dte$variable)
beta <- summary.dte$mean[idx]
cum_haz_ctrl <- design %*% beta

fitted <- gendata
fitted$survival <- NA
fitted$survival[which(fitted$group == "Control")] <- exp(-cum_haz_ctrl)

# Survival curve in treatment arm
# Compute Gaussian quadrature
n_nodes <- 50
gq <- gauss.quad(n = n_nodes, kind = "legendre")
nodes <- gq$nodes
weights <- gq$weights
bdry_knots <- c(0, max(gendata$t))
t1 <- gendata$t[gendata$group == "Treatment"]
design <- iSpline(t1, knots = knots, Boundary.knots = bdry_knots)
df <- ncol(design)
des_gq_1 <- list()
for (l in 1:length(nodes)) {
  input <- t1/2 * (nodes[l] + 1)
  D_l <- iSpline(input, knots=knots, Boundary.knots = bdry_knots, derivs = 1)
  des_gq_1[[l]] <- as.matrix(D_l)
}

tau_mean <- summary.dte %>% filter(variable == "tau") %>% pull(mean)
beta_trt_mean <- summary.dte %>% filter(variable == "beta_trt") %>% pull(mean)
epsilon_mean <- summary.dte %>% filter(variable == "epsilon") %>% pull(mean)

cum_haz_trt <- cumulative_hazard_trt(design_cum = design, t = t1, coef = beta, tau = tau_mean, 
                                     epsilon = epsilon_mean, beta_trt = beta_trt_mean, des_gq = des_gq_1,
                                     nodes = nodes, weights = weights)

fitted$survival[which(fitted$group == "Treatment")] <- exp(-cum_haz_trt)


fit_gendata <- survfit(Surv(time, event) ~ group, data = gendata)
mysurvplot <- ggsurvplot(
  fit_gendata, 
  conf.int = FALSE,
  ggtheme = theme_classic(), 
  data = gendata,
  legend.labs = c("Control", "Treatment"),
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
scale <- 1.8
text_size <- 9*scale
mysurvplot$plot <- mysurvplot$plot + 
  # scale_x_continuous(expand = c(0, 0), breaks = seq(0, 51, 3)) +
  theme(axis.text = element_text(size = 14, color = "black"), 
        legend.text = element_text(size = text_size*0.8, color = "black"), 
        axis.title = element_text(size = text_size),
        axis.line = element_line(color = "black"))

mysurvplot$plot <- mysurvplot$plot + 
  geom_line(data = fitted, aes(time, survival, linetype = group), color = "black", linewidth = 0.85) +
  labs(x = "Time (months)", y = "Probability of Overall Survival", color = "group", linetype = "group") +
  theme(legend.title = element_blank()) +
  scale_linetype_manual(values = c("Control" = "dashed", "Treatment" = "solid"))

text_size <- 9*1.8
scale <- 1.8
mysurvplot$plot <- mysurvplot$plot + 
  # scale_x_continuous(expand = c(0, 0), breaks = seq(0, 51, 3)) +
  theme(axis.text = element_text(size = 14, color = "black"), 
        legend.text = element_text(size = text_size*0.8, color = "black"), 
        axis.title = element_text(size = text_size),
        axis.line = element_line(color = "black"))

pdf("Example/fitted_plot.pdf", height = scale*4, width = scale*7)
mysurvplot$plot / mysurvplot$table + plot_layout(heights = c(8,2))
dev.off()









