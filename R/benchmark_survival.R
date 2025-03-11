#' Benchmark Survival Analysis
#'
#' Performs survival benchmarking across centers using standardized survival ratios and funnel plots.
#'
#' @param data A data frame containing survival data.
#' @param center_col Name of the column with center identifiers.
#' @param status_col Name of the column with event status (0 = censored, 1 = event).
#' @param time_col Name of the column with survival times.
#' @param fixed_time Time point for survival analysis (e.g., 36 for 3 years).
#' @param numeric_covariates Optional vector of column names for numeric covariates used for risk adjustment.
#' @param factor_covariates Optional vector of column names for categorical covariates used for risk adjustment.
#' @param highlight_center Optional center to highlight in the plot.
#' @param conf_levels Confidence levels for funnel plot limits (e.g., c(0.8, 0.95)).
#' @param include_legend Logical; show legend (default TRUE).
#' @param bayesian_adjust Logical; apply Bayesian adjustments (default FALSE).
#' @param bayesian_method Character; "empirical_bayes" or "bayesian_hierarchical" (default "empirical_bayes").
#' @return A ggplot object representing the funnel plot.
#' @importFrom stats as.formula predict qnorm var
#' @importFrom dplyr %>% group_by group_modify summarise ungroup left_join mutate filter across all_of
#' @importFrom ggplot2 ggplot geom_line aes geom_hline geom_point labs scale_color_manual theme_minimal theme
#' @importFrom rlang sym
#' @importFrom tibble tibble
#' @export
#' @examples
#' \dontrun{
#'   data <- data.frame(
#'     center = c("A", "A", "B", "B"),
#'     time = c(10, 20, 15, 25),
#'     status = c(1, 0, 1, 0),
#'     age = c(50, 60, 55, 65)
#'   )
#'   benchmark_survival(data, "center", "status", "time", fixed_time = 15,
#'                      numeric_covariates = "age")
#' }
benchmark_survival <- function(data, center_col, status_col, time_col, fixed_time,
                               numeric_covariates = NULL, factor_covariates = NULL,
                               highlight_center = NULL, conf_levels = c(0.8, 0.95),
                               include_legend = TRUE,
                               bayesian_adjust = FALSE, bayesian_method = "empirical_bayes") {
  # Check dependencies
  if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (bayesian_adjust && bayesian_method == "bayesian_hierarchical" &&
      !requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for bayesian_hierarchical method.")
  }

  # Input validation
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!is.numeric(fixed_time) || fixed_time <= 0) stop("fixed_time must be a positive numeric value")
  if (!bayesian_method %in% c("empirical_bayes", "bayesian_hierarchical")) {
    stop("bayesian_method must be 'empirical_bayes' or 'bayesian_hierarchical'")
  }

  # Process data and fit models
  prepared_data <- prepare_data(data, center_col, status_col, time_col, numeric_covariates, factor_covariates)
  cox_model <- fit_cox_model(prepared_data, numeric_covariates, factor_covariates)
  expected_survival <- compute_expected_survival(prepared_data, cox_model, fixed_time)
  observed_survival <- compute_observed_survival(prepared_data, fixed_time)
  center_summary <- summarize_centers(prepared_data, expected_survival, observed_survival)

  # Bayesian hierarchical adjustment (if applicable)
  if (bayesian_adjust && bayesian_method == "bayesian_hierarchical") {
    bayesian_model <- fit_bayesian_hierarchical(prepared_data)
    metrics <- compute_metrics(center_summary, bayesian_adjust, bayesian_method, bayesian_model)
  } else {
    metrics <- compute_metrics(center_summary, bayesian_adjust, bayesian_method)
  }

  funnel_limits <- generate_funnel_limits(metrics, conf_levels)
  plot <- create_funnel_plot(metrics, funnel_limits, highlight_center, conf_levels, include_legend)

  return(plot)
}

#' Prepare Data for Survival Analysis
#' @keywords internal
prepare_data <- function(data, center_col, status_col, time_col, numeric_covariates, factor_covariates) {
  required_cols <- c(center_col, status_col, time_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

  optional_cols <- c(numeric_covariates, factor_covariates)
  missing_optional <- setdiff(optional_cols, names(data))
  if (length(missing_optional) > 0) {
    warning("Missing optional columns: ", paste(missing_optional, collapse = ", "))
  }

  data <- data %>%
    dplyr::rename(
      center = !!rlang::sym(center_col),
      status = !!rlang::sym(status_col),
      time = !!rlang::sym(time_col)
    ) %>%
    dplyr::mutate(
      status = as.integer(status),
      time = as.numeric(time),
      center = as.factor(center)
    ) %>%
    dplyr::filter(
      !is.na(time),
      !is.na(status),
      time > 0
    )

  if (!is.null(factor_covariates)) {
    data <- data %>% dplyr::mutate(dplyr::across(dplyr::all_of(factor_covariates),
                                                 ~as.factor(ifelse(is.na(.), "Missing", .))))
  }
  if (!is.null(numeric_covariates)) {
    data <- data %>% dplyr::mutate(dplyr::across(dplyr::all_of(numeric_covariates),
                                                 ~ifelse(is.na(.), median(., na.rm = TRUE), .)))
  }

  return(data)
}

#' Fit Cox Proportional Hazards Model
#' @keywords internal
fit_cox_model <- function(data, numeric_covariates, factor_covariates) {
  covariates <- c(numeric_covariates, factor_covariates)
  formula <- if (length(covariates) > 0) {
    stats::as.formula(paste("Surv(time, status) ~", paste(covariates, collapse = " + ")))
  } else {
    stats::as.formula("Surv(time, status) ~ 1")
  }
  survival::coxph(formula, data = data)
}

#' Compute Expected Survival Probabilities
#' @keywords internal
compute_expected_survival <- function(data, cox_model, fixed_time) {
  bh <- survival::basehaz(cox_model, centered = FALSE)
  bh_times <- bh$time
  bh_haz <- bh$hazard
  if (max(bh_times) < fixed_time) {
    warning("fixed_time exceeds maximum observed time; extrapolating hazard.")
    H_fixed <- max(bh_haz)
  } else {
    H_fixed <- stats::approx(bh_times, bh_haz, xout = fixed_time, rule = 2)$y  # Interpolate
  }

  lp <- stats::predict(cox_model, newdata = data, type = "lp")
  data %>%
    dplyr::mutate(S_i = exp(-H_fixed * exp(lp)))
}

#' Compute Observed Survival Probabilities
#' @keywords internal
compute_observed_survival <- function(data, fixed_time) {
  data %>%
    dplyr::group_by(center) %>%
    dplyr::group_modify(~ {
      fit <- survival::survfit(Surv(time, status) ~ 1, data = .x)
      idx <- findInterval(fixed_time, fit$time)
      S_obs <- if (idx == 0 || length(fit$surv) == 0) {
        warning("No events before fixed_time for center ", .y$center, "; setting S_obs to NA.")
        NA_real_
      } else {
        fit$surv[idx]
      }
      tibble::tibble(S_obs = S_obs)
    }) %>%
    dplyr::ungroup()
}

#' Summarize Center-Level Survival Metrics
#' @keywords internal
summarize_centers <- function(data, expected_survival, observed_survival) {
  patient_agg <- expected_survival %>%
    dplyr::group_by(center) %>%
    dplyr::summarise(
      E_S = sum(S_i, na.rm = TRUE),
      V_S = sum(S_i * (1 - S_i), na.rm = TRUE),
      .groups = "drop"
    )

  data %>%
    dplyr::group_by(center) %>%
    dplyr::summarise(n_patients = dplyr::n(), .groups = "drop") %>%
    dplyr::left_join(patient_agg, by = "center") %>%
    dplyr::left_join(observed_survival, by = "center") %>%
    dplyr::mutate(O_S = S_obs * n_patients)
}

#' Fit Bayesian Hierarchical Model
#' @keywords internal
fit_bayesian_hierarchical <- function(data) {
  brms::brm(
    formula = Surv(time, status) ~ 1 + (1 | center),
    data = data,
    family = brms::weibull(),
    prior = c(brms::prior(brms::normal(0, 5), class = "Intercept")),
    chains = 4, iter = 2000, warmup = 500, cores = 4,
    silent = 2  # Suppress verbose output
  )
}

#' Compute Survival Metrics with Optional Bayesian Adjustment
#' @keywords internal
compute_metrics <- function(center_summary, bayesian_adjust, bayesian_method, bayesian_model = NULL) {
  center_summary <- center_summary %>%
    dplyr::mutate(
      SSR = O_S / E_S,
      Z = (O_S - E_S) / sqrt(V_S),
      precision = E_S^2 / V_S
    ) %>%
    dplyr::filter(
      is.finite(SSR),
      is.finite(precision),
      precision > 0
    )

  if (bayesian_adjust) {
    if (bayesian_method == "empirical_bayes") {
      global_mean <- mean(center_summary$SSR, na.rm = TRUE)
      global_var <- stats::var(center_summary$SSR, na.rm = TRUE)
      center_summary <- center_summary %>%
        dplyr::mutate(
          bayes_SSR = (precision * SSR + global_var * global_mean) / (precision + global_var)
        )
    } else if (bayesian_method == "bayesian_hierarchical") {
      if (is.null(bayesian_model)) stop("Bayesian model required for hierarchical method.")
      center_effects <- brms::ranef(bayesian_model)$center[, "Estimate", "Intercept"]
      center_summary <- center_summary %>%
        dplyr::mutate(
          bayes_SSR = SSR * exp(center_effects[as.character(center)])
        )
    }
  } else {
    center_summary <- center_summary %>% dplyr::mutate(bayes_SSR = SSR)
  }

  return(center_summary)
}

#' Generate Funnel Plot Limits
#' @keywords internal
generate_funnel_limits <- function(metrics, conf_levels) {
  precision_seq <- seq(
    max(min(metrics$precision, na.rm = TRUE), 0.01),
    max(metrics$precision, na.rm = TRUE),
    length.out = 100
  )

  limits <- tibble::tibble(precision = precision_seq)
  for (conf in conf_levels) {
    z <- stats::qnorm((1 + conf) / 2)
    limits[[paste0("upper_", conf)]] <- 1 + z / sqrt(precision_seq)
    limits[[paste0("lower_", conf)]] <- pmax(0, 1 - z / sqrt(precision_seq))
  }

  return(limits)
}

#' Create Funnel Plot
#' @keywords internal
create_funnel_plot <- function(metrics, limits, highlight = NULL, conf_levels, include_legend = TRUE) {
  gg <- ggplot2::ggplot() +
    ggplot2::geom_line(data = limits, ggplot2::aes(x = precision, y = upper_0.95), linetype = "dashed", color = "red") +
    ggplot2::geom_line(data = limits, ggplot2::aes(x = precision, y = lower_0.95), linetype = "dashed", color = "red") +
    ggplot2::geom_line(data = limits, ggplot2::aes(x = precision, y = upper_0.8), linetype = "dashed", color = "blue") +
    ggplot2::geom_line(data = limits, ggplot2::aes(x = precision, y = lower_0.8), linetype = "dashed", color = "blue") +
    ggplot2::geom_hline(yintercept = 1, color = "grey50", linetype = 2) +
    ggplot2::geom_point(data = metrics,
                        ggplot2::aes(x = precision, y = bayes_SSR, fill = center,
                                     color = ifelse(center == highlight, "highlighted", "normal")),
                        shape = 21, size = 3) +
    ggplot2::labs(
      x = "Precision (E_S^2 / V_S)",
      y = "Bayesian-Adjusted Standardized Survival Ratio (O/E)",
      title = "Funnel Plot of Survival Ratios by Center"
    ) +
    ggplot2::scale_color_manual(values = c("highlighted" = "red", "normal" = "black"), guide = "none") +
    ggplot2::theme_minimal()

  if (!include_legend) gg <- gg + ggplot2::theme(legend.position = "none")
  return(gg)
}

# Declare global variables to avoid NOTES in R CMD check
utils::globalVariables(c(
  ".", "center", "status", "time", "S_i", "n_patients", "E_S", "V_S", "O_S",
  "SSR", "Z", "precision", "lp", "S_obs", "bayes_SSR", "upper_0.95", "lower_0.95",
  "upper_0.8", "lower_0.8"
))
