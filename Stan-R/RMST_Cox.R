library(survival)
library(dplyr)

# ---------------------------
# Helper: compute RMST from a survfit object up to tau
# sfit: survfit object for a single group/strata
# tau: truncation time
# Returns numeric RMST
# ---------------------------
rmst_from_survfit <- function(sfit, tau) {
  # times we want survival at: 0, (event times < tau), tau
  # Use summary to get survival at necessary times (extend=TRUE to get tail)
  check_times <- sort(unique(c(0, sfit$time[sfit$time < tau], tau)))
  # summary returns a list with $surv; if multiple strata present it returns list per strata
  sumry <- summary(sfit, times = check_times, extend = TRUE)
  # summary(...)$surv will be a vector (no strata) or a matrix if multiple strata;
  survvals <- as.numeric(sumry$surv)
  # integrate using left Riemann (or rectangle) rule: sum_{i} S(t_{i}) * (t_{i+1} - t_{i})
  # but because we included 0, we use survvals[1:(n-1)] * diffs
  times <- check_times
  diffs <- diff(times)
  area <- sum(survvals[1:(length(survvals)-1)] * diffs)
  return(area)
}

# ---- patched helper: ensure newdata row matches model data types ----
.make_template_for_cox <- function(cox_mod, data) {
  # Build a template data.frame matching the model.frame used by coxph.
  mf <- model.frame(delete.response(terms(cox_mod)), data)
  template <- as.data.frame(lapply(mf, function(col) {
    if (is.factor(col)) {
      # set to reference level (first level)
      factor(levels(col)[1], levels = levels(col))
    } else {
      mean(col, na.rm = TRUE)
    }
  }), stringsAsFactors = FALSE)
  # ensure factor columns have same class/levels
  for (nm in names(mf)) {
    if (is.factor(mf[[nm]])) template[[nm]] <- factor(as.character(template[[nm]]), levels = levels(mf[[nm]]))
    else template[[nm]] <- as.numeric(template[[nm]])
  }
  return(template)
}

# ---- patched rmst_diff_cox to be slightly more robust ----
rmst_diff_cox <- function(cox_formula, data, tau) {
  # data must be a plain data.frame here
  data <- as.data.frame(data)
  cox_mod <- coxph(cox_formula, data = data, x = TRUE)
  # build template consistent with model frame
  template <- .make_template_for_cox(cox_mod, data)
  if (!("trt" %in% names(template))) stop("Model/data must contain a variable named 'trt' (0/1).")
  new0 <- template; new0$trt <- 0
  new1 <- template; new1$trt <- 1
  # ensure new0/new1 are data.frames with proper factor levels
  new0 <- as.data.frame(new0, stringsAsFactors = FALSE)
  new1 <- as.data.frame(new1, stringsAsFactors = FALSE)
  for (nm in names(template)) {
    if (is.factor(template[[nm]])) {
      new0[[nm]] <- factor(new0[[nm]], levels = levels(template[[nm]]))
      new1[[nm]] <- factor(new1[[nm]], levels = levels(template[[nm]]))
    }
  }
  
  sfit0 <- survfit(cox_mod, newdata = new0)
  sfit1 <- survfit(cox_mod, newdata = new1)
  
  rmst0 <- rmst_from_survfit(sfit0, tau)
  rmst1 <- rmst_from_survfit(sfit1, tau)
  diff <- rmst1 - rmst0
  
  return(list(cox = cox_mod, rmst0 = rmst0, rmst1 = rmst1, diff = diff))
}

# ---- robust stratified bootstrap + coercion to data.frame ----
bootstrap_rmst_diff_cox <- function(data, cox_formula, tau, nboot = 1000, strata = NULL, seed = 2025, verbose = TRUE) {
  set.seed(seed)
  n <- nrow(data)
  diffs <- numeric(nboot)
  
  for (b in seq_len(nboot)) {
    sampled <- NULL
    # Stratified sampling if requested
    if (!is.null(strata)) {
      if (!(strata %in% names(data))) stop("strata variable not found in data.")
      parts <- split(data, data[[strata]], drop = TRUE)
      sampled_parts <- lapply(parts, function(df_part) {
        df_part[sample.int(nrow(df_part), size = nrow(df_part), replace = TRUE), , drop = FALSE]
      })
      sampled <- do.call(rbind, sampled_parts)
      # keep rownames clean
      rownames(sampled) <- NULL
    } else {
      sampled <- data[sample.int(n, size = n, replace = TRUE), , drop = FALSE]
      rownames(sampled) <- NULL
    }
    
    # Force plain data.frame (no grouped_df, no tibble quirks)
    sampled <- as.data.frame(sampled, stringsAsFactors = FALSE)
    
    # Try compute rmst diff on this bootstrap sample
    diffs[b] <- tryCatch({
      tmp <- rmst_diff_cox(cox_formula, sampled, tau)
      tmp$diff
    }, error = function(e) {
      if (verbose) message(sprintf("bootstrap %d failed: %s", b, e$message))
      NA_real_
    })
  }
  
  valid <- !is.na(diffs)
  if (sum(valid) < nboot * 0.5 && verbose) message(sprintf("Warning: only %d/%d bootstrap replicates succeeded.", sum(valid), nboot))
  diffs_valid <- diffs[valid]
  
  ci_perc <- quantile(diffs_valid, probs = c(0.025, 0.975), na.rm = TRUE)
  se <- sd(diffs_valid, na.rm = TRUE)
  mean_est <- mean(diffs_valid, na.rm = TRUE)
  ci_norm <- mean_est + qnorm(c(0.025, 0.975)) * se
  
  return(list(boot_diffs = diffs_valid,
              nboot_attempted = nboot,
              nboot_success = length(diffs_valid),
              ci_percentile = ci_perc,
              ci_normal = ci_norm,
              se = se))
}





# ## Example Usage
# immuno <- readRDS("/proj/ibrahimlab/dte/Anil/data/ipd_combined_09_25.Rds")
# 
# rmst_est_cox <- rmst_diff_cox(cox_formula = Surv(time, status) ~ trt,
#                               data = immuno %>% mutate(trt = ifelse(trt_group == "Durvalumab + EP", 1, 0)), 
#                               tau = 36)
# 
# rmst_ci_cox <- bootstrap_rmst_diff_cox(data = immuno %>% mutate(trt = ifelse(trt_group == "Durvalumab + EP", 1, 0)), 
#                                        cox_formula = Surv(time, status) ~ trt, tau = 36, nboot = 1000)
# 
# rmst_cox <- c(rmst_est_cox$diff, rmst_ci_cox$ci_percentile, rmst_ci_cox$se)
# names(rmst_cox) <- c("estimate", "lower", "upper", "se")
# 
# ## Also compare to KM method
# library(survRM2)
#  
# rmst_km <- rmst2(immuno$time, immuno$status, ifelse(immuno$trt_group == "Durvalumab + EP", 1, 0), 
#                  tau = 36, alpha = 0.05)
# 
# # This gives difference
# rmst_km <- as.data.frame(t(rmst_km$unadjusted.result[1, 1:3]))
# names(rmst_km) <- c("estimate", "lower", "upper")
# rmst_km$se <- (rmst_km$estimate - rmst_km$lower) / 1.96
# 
# # Save results
# summary <- rbind(rmst_cox, rmst_km)
# summary$method <- c("Cox", "KM")




