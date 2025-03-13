#' Benchmark Survival Analysis
#'
#' Performs survival benchmarking across centers using standardized survival ratios and funnel plots.
#'
#' @param data A data frame containing survival data.
#' @param center_col Name of the column with center identifiers (e.g., "Sitename").
#' @param status_col Name of the column with event status (0 = censored, 1 = event, e.g., "Death").
#' @param time_col Name of the column with survival times (e.g., "DTime").
#' @param fixed_time Time point for survival analysis (e.g., 36 for 36 months).
#' @param highlight_center Optional center to highlight in the plot (e.g., "SiteA").
#' @param conf_levels Confidence levels for funnel plot limits (e.g., c(0.8, 0.95), default c(0.8, 0.95)).
#' @param numeric_covariates Optional vector of column names for numeric covariates (e.g., c("Age_centered")). Set to NULL for no covariates.
#' @param factor_covariates Optional vector of column names for factor covariates (e.g., c("ISS_cat")). Set to NULL for no covariates.
#' @param include_legend Logical; show legend (default TRUE).
#' @param stop_on_ph_violation Logical; stop if PH assumption is violated (default FALSE). Ignored when no covariates are included.
#' @param variance_method Character; method to estimate variance for control limits. Options are "greenwood" (default, uses Greenwood variance from Kaplan-Meier estimator) or "scaling" (uses Greenwood variance with empirical scaling factors).
#' @return A list containing:
#'   \item{plot}{A ggplot object representing the funnel plot.}
#'   \item{metrics}{A data frame with center-level metrics (e.g., SSR, precision, Z).}
#'   \item{funnel_limits}{A data frame with control limits at each precision level.}
#' @export
#' @importFrom dplyr %>% mutate group_by summarise left_join filter rename group_modify ungroup n
#' @importFrom survival coxph basehaz survfit Surv
#' @importFrom ggplot2 ggplot geom_point geom_line labs scale_color_manual theme_minimal coord_cartesian geom_hline scale_linetype_manual aes theme geom_text
#' @importFrom stats as.formula qnorm approx median setNames
#' @importFrom utils globalVariables
#' @importFrom rlang sym .data
#' @importFrom tibble tibble
#' @examples
#' \dontrun{
#'   test_data <- data.frame(
#'     Sitename = rep(c("SiteA", "SiteB"), each = 50),
#'     Death = rbinom(100, 1, 0.3),
#'     DTime = rexp(100, rate = 0.1)
#'   )
#'   result_default <- benchmark_survival(
#'     data = test_data,
#'     center_col = "Sitename",
#'     status_col = "Death",
#'     time_col = "DTime",
#'     fixed_time = 36,
#'     numeric_covariates = NULL,
#'     factor_covariates = NULL,
#'     highlight_center = "SiteA",
#'     conf_levels = c(0.8, 0.95),
#'     include_legend = TRUE
#'   )  # Uses "greenwood" by default
#'   result_scaling <- benchmark_survival(
#'     data = test_data,
#'     center_col = "Sitename",
#'     status_col = "Death",
#'     time_col = "DTime",
#'     fixed_time = 36,
#'     numeric_covariates = NULL,
#'     factor_covariates = NULL,
#'     highlight_center = "SiteA",
#'     conf_levels = c(0.8, 0.95),
#'     include_legend = TRUE,
#'     variance_method = "scaling"
#'   )
#'   print(result_default$plot)
#'   print(result_scaling$plot)
#' }

benchmark_survival <- function(data, center_col, status_col, time_col, fixed_time,
                               highlight_center = NULL, conf_levels = c(0.8, 0.95),
                               numeric_covariates = NULL, factor_covariates = NULL,
                               include_legend = TRUE, stop_on_ph_violation = FALSE,
                               variance_method = c("greenwood", "scaling")) {
  # Input validation
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!is.numeric(fixed_time) || fixed_time <= 0) stop("fixed_time must be a positive numeric value")
  if (!is.null(conf_levels) && !all(conf_levels > 0 & conf_levels < 1)) stop("conf_levels must be between 0 and 1")
  variance_method <- match.arg(variance_method)

  # Processing pipeline
  prepared_data <- prepare_data(data, center_col, status_col, time_col, numeric_covariates, factor_covariates)
  cox_model <- fit_cox_model(prepared_data, c(numeric_covariates, factor_covariates))

  # Proportional Hazards (PH) Assumption Check (skip for null model)
  has_covariates <- !is.null(numeric_covariates) || !is.null(factor_covariates)
  if (has_covariates) {
    ph_test <- survival::cox.zph(cox_model)
    message("\nProportional Hazards (PH) Assumption Test:")
    print(ph_test)

    if (all(ph_test$table[, 3] > 0.05)) {
      message("PH assumption is fulfilled.")
    } else {
      warning("PH assumption may be violated. Consider using time-dependent covariates or stratification.")
      message("Covariates violating PH assumption: ",
              paste(rownames(ph_test$table)[ph_test$table[, 3] <= 0.05], collapse = ", "))
      if (stop_on_ph_violation) stop("PH assumption violated; set stop_on_ph_violation = FALSE to proceed anyway.")
    }
  } else {
    message("\nNo covariates included; skipping Proportional Hazards (PH) Assumption Test.")
  }

  expected_survival <- compute_expected_survival(prepared_data, cox_model, fixed_time, c(numeric_covariates, factor_covariates))
  observed_survival <- compute_observed_survival(prepared_data, fixed_time, variance_method)
  center_summary <- summarize_centers(prepared_data, expected_survival, observed_survival, variance_method)
  metrics <- compute_metrics(center_summary, variance_method)
  funnel_limits <- generate_funnel_limits(metrics, conf_levels, variance_method)

  if (is.null(funnel_limits)) {
    stop("Unable to generate funnel plot limits due to invalid precision values.")
  }

  plot <- create_funnel_plot(metrics, funnel_limits, highlight_center, conf_levels, include_legend)

  # Model Fit Assessment
  cox_fit_summary <- summary(cox_model)
  message("\nModel Fit Assessment:")
  message("Concordance Index (C-statistic):")
  print(cox_fit_summary$concordance)
  message("Likelihood Ratio Test:")
  print(cox_fit_summary$logtest)
  message("Wald Test:")
  print(cox_fit_summary$waldtest)
  message("Score (Log-Rank) Test:")
  print(cox_fit_summary$sctest)

  # Return a list with plot, metrics, and funnel_limits
  return(list(
    plot = plot,
    metrics = metrics,
    funnel_limits = funnel_limits
  ))
}

#' @rdname benchmark_survival
#' @keywords internal
prepare_data <- function(data, center_col, status_col, time_col, numeric_covariates, factor_covariates) {
  required_cols <- c(center_col, status_col, time_col)
  if (!is.null(numeric_covariates)) required_cols <- c(required_cols, numeric_covariates)
  if (!is.null(factor_covariates)) required_cols <- c(required_cols, factor_covariates)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

  data <- data %>%
    rename(
      center = !!sym(center_col),
      status = !!sym(status_col),
      time   = !!sym(time_col)
    ) %>%
    mutate(
      center = as.factor(center),
      status = as.integer(status),
      time = as.numeric(time)
    ) %>%
    filter(
      !is.na(time),
      !is.na(status),
      time > 0
    )

  # Handle numeric covariates
  if (!is.null(numeric_covariates)) {
    for (cov in numeric_covariates) {
      if (!is.numeric(data[[cov]])) {
        warning("Converting numeric covariate '", cov, "' to numeric.")
        data[[cov]] <- as.numeric(as.character(data[[cov]]))
      }
      # Impute missing values with median
      if (any(is.na(data[[cov]]))) {
        data[[cov]][is.na(data[[cov]])] <- median(data[[cov]], na.rm = TRUE)
      }
    }
  }

  # Handle factor covariates
  if (!is.null(factor_covariates)) {
    for (cov in factor_covariates) {
      if (!is.factor(data[[cov]])) {
        warning("Converting factor covariate '", cov, "' to factor.")
        data[[cov]] <- as.factor(data[[cov]])
      }
      # Impute missing values with a "Missing" category
      if (any(is.na(data[[cov]]))) {
        data[[cov]][is.na(data[[cov]])] <- "Missing"
        data[[cov]] <- as.factor(data[[cov]])
      }
      # Ensure all levels are present in the data
      data[[cov]] <- droplevels(data[[cov]])
    }
  }

  return(data)
}

#' @rdname benchmark_survival
#' @keywords internal
fit_cox_model <- function(data, covariates) {
  if (is.null(covariates) || length(covariates) == 0) {
    formula <- as.formula("Surv(time, status) ~ 1")
  } else {
    formula <- as.formula(paste("Surv(time, status) ~", paste(covariates, collapse = " + ")))
  }
  survival::coxph(formula, data = data)
}

#' @rdname benchmark_survival
#' @keywords internal
compute_expected_survival <- function(data, cox_model, fixed_time, covariates) {
  bh <- survival::basehaz(cox_model, centered = FALSE)
  bh_times <- bh$time
  bh_haz <- bh$hazard

  # Interpolate or extrapolate the baseline hazard at fixed_time
  if (max(bh_times) < fixed_time) {
    warning("fixed_time exceeds maximum observed time; extrapolating hazard.")
    H_fixed <- max(bh_haz)
  } else {
    H_fixed <- stats::approx(bh_times, bh_haz, xout = fixed_time, rule = 2)$y
  }

  # Compute linear predictor using predict
  if (is.null(covariates) || length(covariates) == 0) {
    lp <- rep(0, nrow(data))
  } else {
    lp <- predict(cox_model, newdata = data, type = "lp")
  }

  data %>%
    mutate(
      S_i = exp(-H_fixed * exp(lp))
    )
}

#' @rdname benchmark_survival
#' @keywords internal
compute_observed_survival <- function(data, fixed_time, variance_method) {
  data %>%
    group_by(center) %>%
    group_modify(~ {
      fit <- survival::survfit(Surv(time, status) ~ 1, data = .x)
      idx <- findInterval(fixed_time, fit$time)
      S_obs <- if (idx == 0 || length(fit$surv) == 0) {
        warning("No events before fixed_time for center ", .y$center, "; setting S_obs to 1.")
        1
      } else {
        fit$surv[idx]
      }
      # Greenwood variance for S_obs
      var_S_obs <- if (idx == 0 || length(fit$std.err) == 0) {
        0  # No variance if no data
      } else {
        (fit$std.err[idx] * S_obs)^2  # Variance of S_obs using Greenwood formula
      }
      tibble(S_obs = S_obs, var_S_obs = var_S_obs)
    }) %>%
    ungroup()
}

#' @rdname benchmark_survival
#' @keywords internal
summarize_centers <- function(data, expected_survival, observed_survival, variance_method) {
  patient_agg <- expected_survival %>%
    group_by(center) %>%
    summarise(
      E_S = sum(S_i, na.rm = TRUE),
      V_S = sum(S_i * (1 - S_i), na.rm = TRUE),
      .groups = "drop"
    )

  data %>%
    group_by(center) %>%
    summarise(n_patients = n(), .groups = "drop") %>%
    left_join(patient_agg, by = "center") %>%
    left_join(observed_survival, by = "center") %>%
    mutate(
      O_S = S_obs * n_patients,
      V_O_S = var_S_obs * (n_patients^2)  # Always use Greenwood variance for O_S
    )
}

#' @rdname benchmark_survival
#' @keywords internal
compute_metrics <- function(center_summary, variance_method) {
  center_summary %>%
    mutate(
      SSR = O_S / E_S,
      Z = (O_S - E_S) / sqrt(V_O_S),
      precision = E_S^2 / V_O_S
    ) %>%
    filter(
      is.finite(SSR),
      is.finite(precision),
      precision > 0,
      !is.na(E_S) & E_S > 0,
      !is.na(V_O_S) & V_O_S > 0
    )
}

#' @rdname benchmark_survival
#' @keywords internal
generate_funnel_limits <- function(metrics, conf_levels, variance_method) {
  # Check for valid precision values
  valid_precision <- metrics$precision[!is.na(metrics$precision) & is.finite(metrics$precision)]
  if (length(valid_precision) < 2) {
    warning("Not enough valid precision values to generate funnel limits.")
    return(NULL)
  }

  precision_seq <- seq(
    max(min(valid_precision, na.rm = TRUE), 0.01),
    max(valid_precision, na.rm = TRUE),
    length.out = 100
  )

  limits <- tibble(precision = precision_seq)

  if (variance_method == "scaling") {
    # Apply scaling factors to adjust control limits
    scaling_80 <- 0.654  # Adjusted for 0.0503 observed vs 0.2 expected
    scaling_95 <- 0.781  # Adjusted for 0.0121 observed vs 0.05 expected
  } else {
    scaling_80 <- 1
    scaling_95 <- 1
  }

  for (i in seq_along(conf_levels)) {
    conf <- conf_levels[i]
    z <- qnorm((1 + conf) / 2)
    scaling <- if (conf == 0.8) scaling_80 else if (conf == 0.95) scaling_95 else 1
    z_adjusted <- z * scaling
    limits[[paste0("upper_", conf)]] <- 1 + z_adjusted / sqrt(precision_seq)
    limits[[paste0("lower_", conf)]] <- pmax(0, 1 - z_adjusted / sqrt(precision_seq))
  }

  limits
}

#' @rdname benchmark_survival
#' @keywords internal
create_funnel_plot <- function(metrics, limits, highlight = NULL, conf_levels, include_legend = TRUE) {
  # Define color and linetype mappings for confidence levels
  color_map <- stats::setNames(c("blue", "red"), c("80%", "95%"))
  linetype_map <- stats::setNames(c("dashed", "dotted"), c("80%", "95%"))
  y_min <- min(metrics$SSR, na.rm = TRUE) - 0.1  # Adjust y-axis limit

  gg <- ggplot() +
    # Plot centers as points with fill = center (shape = 21 allows separate fill and outline)
    geom_point(
      data = metrics,
      aes(x = precision, y = SSR, fill = center),
      shape = 21,
      color = "black",
      size = 3,
      alpha = 0.7
    ) +
    # Horizontal reference line at SSR=1
    geom_hline(yintercept = 1, color = "grey50", linetype = 2) +
    theme_minimal(base_size = 14) +
    labs(
      title  = "Survival Outcomes Benchmarking",
      x      = "Precision (E^2/V)",
      y      = "Standardized Survival Ratio (O/E)",
      fill   = "Center",
      color  = "Confidence",
      linetype = "Confidence"
    ) +
    scale_color_manual(values = color_map, drop = FALSE) +
    scale_linetype_manual(values = linetype_map, drop = FALSE) +
    coord_cartesian(ylim = c(y_min, NA)) +
    theme(legend.position = if (include_legend) "right" else "none",
          legend.justification = c(1, 1),
          legend.box = "vertical")

  # Plot upper and lower funnel lines for each confidence level
  for (conf in conf_levels) {
    conf_label <- paste0(conf * 100, "%")

    upper_df <- tibble(
      precision = limits$precision,
      SSR       = limits[[paste0("upper_", conf)]],
      conf      = conf_label
    )
    lower_df <- tibble(
      precision = limits$precision,
      SSR       = limits[[paste0("lower_", conf)]],
      conf      = conf_label
    )

    gg <- gg +
      geom_line(
        data = upper_df,
        aes(x = precision, y = SSR, color = conf, linetype = conf),
        size = 1
      ) +
      geom_line(
        data = lower_df,
        aes(x = precision, y = SSR, color = conf, linetype = conf),
        size = 1
      )
  }

  # Highlight and label a center if specified
  if (!is.null(highlight)) {
    highlighted_data <- dplyr::filter(metrics, center == highlight)
    if (nrow(highlighted_data) > 0) {
      gg <- gg +
        geom_point(
          data = highlighted_data,
          aes(x = precision, y = SSR),
          shape = 21,
          fill = "red",
          color = "black",
          size = 5,
          stroke = 1.3,
          show.legend = FALSE
        ) +
        geom_text(
          data = highlighted_data,
          aes(x = precision, y = SSR, label = center),
          vjust = -1,
          color = "red",
          size = 5,
          fontface = "bold"
        )
    }
  }

  return(gg)
}

utils::globalVariables(c(
  ".", "center", "status", "time", "S_i", "n_patients",
  "E_S", "V_S", "O_S", "SSR", "Z", "precision", "S_obs",
  "conf", "predict", "var_S_obs", "V_O_S"
))
