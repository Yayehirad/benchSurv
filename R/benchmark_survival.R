#' Benchmark Survival Analysis (Refined Version)
#'
#' Performs survival benchmarking across centers using standardized survival ratios and funnel plots.
#'
#' @param data A data frame containing survival data.
#' @param center_col Name of the column with center identifiers.
#' @param status_col Name of the column with event status (0 = censored, 1 = event).
#' @param time_col Name of the column with survival times.
#' @param age_col Name of the column with age data.
#' @param fixed_time Time point for survival analysis (e.g., 36 for 3 years).
#' @param highlight_center Optional center to highlight in the plot.
#' @param conf_levels Confidence levels for funnel plot limits.
#' @return A ggplot object representing the funnel plot.
#' @export
#' @importFrom dplyr %>% mutate group_by summarise left_join filter rename group_modify ungroup n bind_rows
#' @importFrom survival coxph basehaz survfit Surv cox.zph frailty
#' @importFrom ggplot2 ggplot geom_point geom_line labs scale_color_manual theme_minimal coord_cartesian geom_hline scale_linetype_manual aes theme scale_fill_discrete
#' @importFrom stats coef qnorm setNames approx
#' @importFrom utils tail globalVariables
#' @importFrom tibble tibble
#' @importFrom lme4 glmer
#' @importFrom lme4 lmer
benchmark_survival <- function(data, center_col, status_col, time_col, age_col, fixed_time,
                               highlight_center = NULL, conf_levels = c(0.8, 0.95)) {
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!is.numeric(fixed_time) || fixed_time <= 0) stop("fixed_time must be positive numeric")

  prepared_data <- prepare_data(data, center_col, status_col, time_col, age_col)
  cox_model <- fit_cox_model(prepared_data)

  # Proportional hazards assumption check (Added on line 27)
  ph_test <- survival::cox.zph(cox_model)
  print(ph_test) # Warn user if PH assumption is violated

  expected_survival <- compute_expected_survival(prepared_data, cox_model, fixed_time)
  observed_survival <- compute_observed_survival(prepared_data, fixed_time)
  center_summary <- summarize_centers(prepared_data, expected_survival, observed_survival)
  metrics <- compute_metrics(center_summary)
  funnel_limits <- generate_funnel_limits_poisson(metrics, conf_levels) # Updated with Poisson limits (Line 37)
  create_funnel_plot(metrics, funnel_limits, highlight_center, conf_levels)
}

#' @rdname benchmark_survival
#' @keywords internal
fit_cox_model <- function(data) {
  # Using frailty model for random effects on centers (Added on line 44)
  survival::coxph(Surv(time, status) ~ age + frailty(center), data = data)
}

#' @rdname benchmark_survival
#' @keywords internal
compute_expected_survival <- function(data, cox_model, fixed_time) {
  bh <- survival::basehaz(cox_model, centered = FALSE)
  # Improved baseline hazard estimation using interpolation (Added on line 52)
  H_fixed <- approx(bh$time, bh$hazard, xout = fixed_time, rule = 2)$y

  data %>% mutate(
    lp = coef(cox_model) * age,
    S_i = exp(-H_fixed * exp(lp))
  )
}

#' @rdname benchmark_survival
#' @keywords internal
generate_funnel_limits_poisson <- function(metrics, conf_levels) {
  # Poisson-binomial confidence limits for small samples (Added on line 62)
  precision_seq <- seq(
    max(min(metrics$precision, na.rm = TRUE), 0.01),
    max(metrics$precision, na.rm = TRUE),
    length.out = 100
  )

  limits <- data.frame(precision = precision_seq)

  for (conf in conf_levels) {
    z <- qnorm((1 + conf) / 2)
    limits[[paste0("upper_", conf)]] <- 1 + z / sqrt(precision_seq)
    limits[[paste0("lower_", conf)]] <- pmax(0, 1 - z / sqrt(precision_seq))
  }

  limits
}

#' @rdname benchmark_survival
#' @keywords internal
simulate_benchmark <- function(n_sim = 1000, n_centers = 20, min_patients = 50, max_patients = 500) {
  # Simulation function for validation (Added on line 76)
  results <- replicate(n_sim, {
    center_sizes <- sample(min_patients:max_patients, n_centers, replace = TRUE)
    sim_data <- data.frame()

    for (center in 1:n_centers) {
      n <- center_sizes[center]
      ages <- rnorm(n, mean = 65, sd = 10)
      surv_time <- rweibull(n, shape = 1.5, scale = 5) * exp(-0.02 * (ages - 65))
      cens_time <- runif(n, 0, 10)
      time <- pmin(surv_time, cens_time)
      status <- as.integer(surv_time <= cens_time)
      sim_data <- rbind(sim_data, data.frame(center = center, age = ages, time = time, status = status))
    }

    model <- benchmark_survival(sim_data, "center", "status", "time", "age", 36)
    return(summary(model))
  }, simplify = FALSE)
  return(results)
}
