#' Benchmark Survival Analysis
#'
#' Performs survival benchmarking across centers using standardized survival ratios and funnel plots.
#'
#' @param data A data frame containing survival data.
#' @param center_col Name of the column with center identifiers.
#' @param status_col Name of the column with event status (0 = censored, 1 = event).
#' @param time_col Name of the column with survival times.
#' @param fixed_time Time point for survival analysis (e.g., 36 for 3 years).
#' @param covariates Optional vector of column names for covariates used for risk adjustment.
#' @param numeric_covariates Optional vector of column names for numeric covariates used for risk adjustment.
#' @param factor_covariates Optional vector of column names for categorical covariates used for risk adjustment.
#' @param highlight_center Optional center to highlight in the plot.
#' @param conf_levels Confidence levels for funnel plot limits.
#' @param include_legend Logical; show legend (default TRUE).
#' @param bayesian_adjust Logical; apply Bayesian adjustments (default FALSE).
#' @param bayesian_method Character; "empirical_bayes" or "bayesian_hierarchical" (default "empirical_bayes").
#' @return A ggplot object representing the funnel plot.
#' @export
benchmark_survival <- function(data, center_col, status_col, time_col, fixed_time,
                               covariates = NULL, numeric_covariates = NULL, factor_covariates = NULL,
                               highlight_center = NULL, conf_levels = c(0.8, 0.95),
                               include_legend = TRUE,
                               bayesian_adjust = FALSE, bayesian_method = "empirical_bayes") {
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!is.numeric(fixed_time) || fixed_time <= 0) stop("fixed_time must be positive numeric")

  prepared_data <- prepare_data(data, center_col, status_col, time_col, numeric_covariates, factor_covariates)
  cox_model <- fit_cox_model(prepared_data, numeric_covariates, factor_covariates)

  expected_survival <- compute_expected_survival(prepared_data, cox_model, fixed_time)
  observed_survival <- compute_observed_survival(prepared_data, fixed_time)
  center_summary <- summarize_centers(prepared_data, expected_survival, observed_survival)
  metrics <- compute_metrics(center_summary, bayesian_adjust, bayesian_method)
  funnel_limits <- generate_funnel_limits(metrics, conf_levels)
  plot <- create_funnel_plot(metrics, funnel_limits, highlight_center, conf_levels, include_legend)

  return(plot)
}

#' @rdname benchmark_survival
#' @keywords internal
prepare_data <- function(data, center_col, status_col, time_col, numeric_covariates, factor_covariates) {
  required_cols <- c(center_col, status_col, time_col, numeric_covariates, factor_covariates)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

  data <- data %>%
    rename(
      center = !!sym(center_col),
      status = !!sym(status_col),
      time   = !!sym(time_col)
    ) %>%
    mutate(
      status = as.integer(status),
      time = as.numeric(time)
    ) %>%
    filter(
      !is.na(time),
      !is.na(status),
      time > 0
    )

  if (!is.null(factor_covariates)) {
    data <- data %>% mutate(across(all_of(factor_covariates), ~as.factor(ifelse(is.na(.), "Missing", .))))
  }

  return(data)
}

#' @rdname benchmark_survival
#' @keywords internal
compute_metrics <- function(center_summary, bayesian_adjust, bayesian_method) {
  center_summary <- center_summary %>%
    mutate(
      SSR = O_S / E_S,
      Z = (O_S - E_S) / sqrt(V_S),
      precision = E_S^2 / V_S
    ) %>%
    filter(
      is.finite(SSR),
      is.finite(precision),
      precision > 0
    )

  if (bayesian_adjust) {
    if (bayesian_method == "empirical_bayes") {
      global_mean <- mean(center_summary$SSR, na.rm = TRUE)
      global_var  <- var(center_summary$SSR, na.rm = TRUE)
      center_summary <- center_summary %>%
        mutate(
          bayes_SSR = (precision * SSR + global_var * global_mean) / (precision + global_var)
        )
    } else if (bayesian_method == "bayesian_hierarchical") {
      library(brms)
      bayesian_model <- brm(
        formula = Surv(time, status) ~ (1 | center),
        data = center_summary,
        family = weibull(),
        prior = c(prior(normal(0, 5), class = "Intercept")),
        chains = 4, iter = 2000, warmup = 500, cores = 4
      )
      center_summary <- center_summary %>%
        mutate(bayes_SSR = fitted(bayesian_model, newdata = center_summary)[, "Estimate"])
    }
  }

  return(center_summary)
}

#' @rdname benchmark_survival
#' @keywords internal
generate_funnel_limits <- function(metrics, conf_levels) {
  precision_seq <- seq(
    max(min(metrics$precision, na.rm = TRUE), 0.01),
    max(metrics$precision, na.rm = TRUE),
    length.out = 100
  )

  limits <- tibble(precision = precision_seq)

  for (conf in conf_levels) {
    z <- qnorm((1 + conf) / 2)
    limits[[paste0("upper_", conf)]] <- 1 + z / sqrt(precision_seq)
    limits[[paste0("lower_", conf)]] <- pmax(0, 1 - z / sqrt(precision_seq))
  }

  limits
}

#' @rdname benchmark_survival
#' @keywords internal
create_funnel_plot <- function(metrics, limits, highlight = NULL, conf_levels, include_legend = TRUE) {
  gg <- ggplot() +
    geom_point(data = metrics, aes(x = precision, y = bayes_SSR, fill = center), shape = 21, color = "black") +
    geom_hline(yintercept = 1, color = "grey50", linetype = 2) +
    labs(y = "Bayesian-Adjusted Standardized Survival Ratio (O/E)")

  return(gg)
}

utils::globalVariables(c(
  ".", "center", "status", "time", "age", "S_i", "n_patients",
  "E_S", "V_S", "O_S", "SSR", "Z", "precision", "lp", "S_obs"
))
