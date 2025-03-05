#' Benchmark Survival Analysis
#'
#' Performs survival benchmarking across centers using standardized survival ratios and funnel plots.
#'
#' @param data A data frame containing survival data.
#' @param center_col Name of the column with center identifiers.
#' @param status_col Name of the column with event status (0 = censored, 1 = event).
#' @param time_col Name of the column with survival times.
#' @param fixed_time Time point for survival analysis (e.g., 36 for 3 years).
#' @param numeric_covariates Vector of column names for numeric covariates used for risk adjustment.
#' @param factor_covariates Vector of column names for categorical covariates used for risk adjustment.
#' @param time_dep_covariates Optional vector of numeric covariate names that violate the PH assumption.
#'   For each, an interaction with log(fixed_time) is included in the model.
#' @param strata_covariates Optional vector of column names to be used for stratification via strata().
#' @param highlight_center Optional center to highlight in the plot.
#' @param conf_levels Confidence levels for funnel plot limits.
#' @param include_legend Logical; show legend (default TRUE).
#' @return A ggplot object representing the funnel plot.
#' @export
#' @importFrom dplyr %>% mutate group_by summarise left_join filter rename group_modify ungroup n bind_rows
#' @importFrom survival coxph basehaz survfit Surv cox.zph
#' @importFrom ggplot2 ggplot geom_point geom_line labs scale_color_manual theme_minimal coord_cartesian geom_hline scale_linetype_manual aes theme scale_fill_discrete
#' @importFrom colorspace rainbow_hcl
#' @importFrom stats coef qnorm setNames
#' @importFrom utils tail globalVariables
#' @importFrom rlang sym .data
#' @importFrom tibble tibble
benchmark_survival <- function(data, center_col, status_col, time_col, fixed_time,
                               numeric_covariates = NULL, factor_covariates = NULL,
                               time_dep_covariates = NULL, strata_covariates = NULL,
                               highlight_center = NULL, conf_levels = c(0.8, 0.95),
                               include_legend = TRUE) {
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!is.numeric(fixed_time) || fixed_time <= 0) stop("fixed_time must be positive numeric")

  # Combine risk adjustment covariates that are not modeled with time-dependency
  base_covariates <- c(numeric_covariates, factor_covariates)

  # Processing pipeline: include strata_covariates in data preparation so we ensure they exist
  prepared_data <- prepare_data(data, center_col, status_col, time_col,
                                covariates = base_covariates,
                                factor_covariates = factor_covariates,
                                strata_covariates = strata_covariates)

  cox_model <- fit_cox_model(prepared_data, base_covariates, strata_covariates,
                             time_dep_covariates, fixed_time)

  # Proportional Hazards (PH) Assumption Check
  ph_test <- survival::cox.zph(cox_model)
  message("\nProportional Hazards (PH) Assumption Test:")
  print(ph_test)

  if (all(ph_test$table[, 3] > 0.05)) {
    message("✅ The Proportional Hazards (PH) assumption is fulfilled.")
  } else {
    warning("⚠️ PH assumption may be violated. Consider modeling the violating covariates as time-dependent.")
    message("Covariates violating PH assumption: ",
            paste(rownames(ph_test$table)[ph_test$table[, 3] <= 0.05], collapse = ", "))
  }

  expected_survival <- compute_expected_survival(prepared_data, cox_model, fixed_time)
  observed_survival <- compute_observed_survival(prepared_data, fixed_time)
  center_summary <- summarize_centers(prepared_data, expected_survival, observed_survival)
  metrics <- compute_metrics(center_summary)
  funnel_limits <- generate_funnel_limits(metrics, conf_levels)
  plot <- create_funnel_plot(metrics, funnel_limits, highlight_center, conf_levels, include_legend)

  # Model Fit Assessment
  cox_fit_summary <- summary(cox_model)
  message("\nModel Fit Assessment:")
  message("✔️ Concordance Index (C-statistic):")
  print(cox_fit_summary$concordance)
  message("✔️ Likelihood Ratio Test:")
  print(cox_fit_summary$logtest)
  message("✔️ Wald Test:")
  print(cox_fit_summary$waldtest)
  message("✔️ Score (Log-Rank) Test:")
  print(cox_fit_summary$sctest)

  return(plot)
}

#' @rdname benchmark_survival
#' @keywords internal
prepare_data <- function(data, center_col, status_col, time_col,
                         covariates, factor_covariates, strata_covariates = NULL) {
  # Ensure required columns exist: add strata_covariates if provided
  required_cols <- unique(c(center_col, status_col, time_col, covariates, strata_covariates))
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

  # Convert specified factor covariates to factors
  if (!is.null(factor_covariates)) {
    for (var in factor_covariates) {
      if (var %in% names(data)) {
        data[[var]] <- as.factor(data[[var]])
      }
    }
  }

  # Also convert strata_covariates to factors if not numeric
  if (!is.null(strata_covariates)) {
    for (var in strata_covariates) {
      if (var %in% names(data) && !is.numeric(data[[var]])) {
        data[[var]] <- as.factor(data[[var]])
      }
    }
  }

  return(data)
}

#' @rdname benchmark_survival
#' @keywords internal
fit_cox_model <- function(data, covariates, strata_covariates = NULL,
                          time_dep_covariates = NULL, fixed_time) {
  # Build the risk adjustment part of the formula from base (non-time-dependent) covariates.
  risk_terms <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else "1"

  # Add time-dependent covariates as both the baseline and an interaction with log(fixed_time)
  if (!is.null(time_dep_covariates)) {
    td_terms <- paste0(time_dep_covariates, " + I(", time_dep_covariates, " * log(", fixed_time, "))",
                       collapse = " + ")
    risk_terms <- paste(risk_terms, td_terms, sep = " + ")
  }

  # Append strata() terms if provided
  if (!is.null(strata_covariates)) {
    strata_formula <- paste0("strata(", paste(strata_covariates, collapse = "), strata("), ")")
    full_formula <- as.formula(paste("Surv(time, status) ~", risk_terms, "+", strata_formula))
  } else {
    full_formula <- as.formula(paste("Surv(time, status) ~", risk_terms))
  }

  survival::coxph(full_formula, data = data)
}

#' @rdname benchmark_survival
#' @keywords internal
compute_expected_survival <- function(data, cox_model, fixed_time) {
  bh <- survival::basehaz(cox_model, centered = FALSE)
  H_fixed <- max(bh$hazard[bh$time <= fixed_time], 0)

  # Using predict() to compute the linear predictor which now includes any time-dependent terms.
  lp <- predict(cox_model, newdata = data, type = "lp")

  data %>%
    mutate(
      S_i = exp(-H_fixed * exp(lp))
    )
}

compute_observed_survival <- function(data, fixed_time) {
  data %>%
    group_by(center) %>%
    group_modify(~ {
      fit <- survival::survfit(Surv(time, status) ~ 1, data = .x)
      idx <- findInterval(fixed_time, fit$time)
      S_obs <- if (idx == 0) 1 else fit$surv[idx]
      tibble(S_obs = S_obs)
    }) %>%
    ungroup()
}

summarize_centers <- function(data, expected_survival, observed_survival) {
  patient_agg <- expected_survival %>%
    group_by(center) %>%
    summarise(
      E_S = sum(S_i),
      V_S = sum(S_i * (1 - S_i)),
      .groups = "drop"
    )

  data %>%
    group_by(center) %>%
    summarise(n_patients = n(), .groups = "drop") %>%
    left_join(patient_agg, by = "center") %>%
    left_join(observed_survival, by = "center") %>%
    mutate(O_S = S_obs * n_patients)
}

compute_metrics <- function(center_summary) {
  center_summary %>%
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
}

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

create_funnel_plot <- function(metrics, limits, highlight = NULL, conf_levels, include_legend = TRUE) {
  # Define color and linetype mappings for confidence levels
  color_map <- c("80%" = "blue", "95%" = "red")
  linetype_map <- c("80%" = "dashed", "95%" = "dotted")
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
    ggplot2::scale_fill_discrete() +
    ggplot2::scale_color_manual(values = color_map, drop = FALSE) +
    scale_linetype_manual(values = linetype_map, drop = FALSE) +
    coord_cartesian(ylim = c(y_min, NA)) +
    theme(legend.position = if (include_legend) "right" else "none",
          legend.justification = c(1, 1),
          legend.box = "vertical")

  # Plot upper and lower funnel lines for each confidence level
  for (conf in conf_levels) {
    conf_label <- paste0(conf * 100, "%")

    upper_df <- tibble::tibble(
      precision = limits$precision,
      SSR       = limits[[paste0("upper_", conf)]],
      conf      = conf_label
    )
    lower_df <- tibble::tibble(
      precision = limits$precision,
      SSR       = limits[[paste0("lower_", conf)]],
      conf      = conf_label
    )

    gg <- gg + geom_line(
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
      ggplot2::geom_text(
        data = highlighted_data,
        aes(x = precision, y = SSR, label = center),
        vjust = -1,    # position label above the point
        color = "red",
        size = 5,
        fontface = "bold"
      )
  }

  return(gg)
}

utils::globalVariables(c(
  ".", "center", "status", "time", "S_i", "n_patients",
  "E_S", "V_S", "O_S", "SSR", "Z", "precision", "lp", "S_obs"
))
