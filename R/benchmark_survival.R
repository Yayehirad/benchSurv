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
#'   If provided and neither numeric_covariates nor factor_covariates are specified,
#'   these covariates will be treated as categorical variables.
#' @param numeric_covariates Optional vector of column names for numeric covariates used for risk adjustment (e.g., c("Age")).
#' @param factor_covariates Optional vector of column names for categorical covariates used for risk adjustment (e.g., c("ISS")).
#' @param highlight_center Optional center to highlight in the plot.
#' @param conf_levels Confidence levels for funnel plot limits.
#' @param include_legend Logical; show legend (default TRUE).
#' @return A ggplot object representing the funnel plot.
#' @export
#' @importFrom dplyr %>% mutate group_by summarise left_join filter rename group_modify ungroup n bind_rows
#' @importFrom survival coxph basehaz survfit Surv cox.zph
#' @importFrom ggplot2 ggplot geom_point geom_line labs scale_color_manual theme_minimal coord_cartesian geom_hline scale_linetype_manual aes theme scale_fill_discrete geom_text
#' @importFrom colorspace rainbow_hcl
#' @importFrom stats coef qnorm setNames
#' @importFrom utils tail globalVariables
#' @importFrom rlang sym .data
#' @importFrom tibble tibble
benchmark_survival <- function(data, center_col, status_col, time_col, fixed_time,
                               covariates = NULL,                 # MODIFIED: added covariates parameter
                               numeric_covariates = NULL, factor_covariates = NULL,
                               highlight_center = NULL, conf_levels = c(0.8, 0.95),
                               include_legend = TRUE) {
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!is.numeric(fixed_time) || fixed_time <= 0) stop("fixed_time must be positive numeric")

  # MODIFIED: If covariates is provided (and numeric/factor not specified), use them as factor_covariates
  if (!is.null(covariates)) {
    if (is.null(numeric_covariates) && is.null(factor_covariates)) {
      factor_covariates <- covariates
    } else {
      warning("Both covariates and numeric/factor covariates provided, ignoring covariates.")
    }
  }

  # Processing pipeline
  prepared_data <- prepare_data(data, center_col, status_col, time_col, numeric_covariates, factor_covariates)
  cox_model <- fit_cox_model(prepared_data, numeric_covariates, factor_covariates)

  # Proportional Hazards (PH) Assumption Check
  if (length(coef(cox_model)) > 0) {   # MODIFIED: Check if model is non-null
    ph_test <- survival::cox.zph(cox_model)
    message("\nProportional Hazards (PH) Assumption Test:")
    print(ph_test)
    if (all(ph_test$table[, 3] > 0.05)) {
      message("✅ The Proportional Hazards (PH) assumption is fulfilled.")
    } else {
      warning("⚠️ PH assumption may be violated. Consider using time-dependent covariates or stratification.")
      message("Covariates violating PH assumption: ",
              paste(rownames(ph_test$table)[ph_test$table[, 3] <= 0.05], collapse = ", "))
    }
  } else {
    message("\nPH test skipped: Null model (no covariates).")  # MODIFIED: Skip PH test if model is null
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
    # Convert factor covariates to factors and replace NA with "Missing" level
    data <- data %>% mutate(across(all_of(factor_covariates), ~as.factor(ifelse(is.na(.), "Missing", .))))
  }

  return(data)
}

#' @rdname benchmark_survival
#' @keywords internal
fit_cox_model <- function(data, numeric_covariates, factor_covariates) {
  covariates <- c(numeric_covariates, factor_covariates)
  formula <- if (length(covariates) > 0) {
    as.formula(paste("Surv(time, status) ~", paste(covariates, collapse = " + ")))
  } else {
    as.formula("Surv(time, status) ~ 1")
  }
  survival::coxph(formula, data = data)
}

#' @rdname benchmark_survival
#' @keywords internal
# MODIFIED: Removed extra covariate parameters and use predict() to compute the linear predictor
compute_expected_survival <- function(data, cox_model, fixed_time) {
  bh <- survival::basehaz(cox_model, centered = FALSE)
  H_fixed <- max(bh$hazard[bh$time <= fixed_time], 0)

  lp <- predict(cox_model, newdata = data, type = "lp")  # MODIFIED: Use predict() for linear predictor

  data %>%
    mutate(
      S_i = exp(-H_fixed * exp(lp))
    )
}

#' @rdname benchmark_survival
#' @keywords internal
compute_observed_survival <- function(data, fixed_time) {
  data %>%
    group_by(center) %>%
    group_modify(~ {
      fit <- survival::survfit(Surv(time, status) ~ 1, data = .x)
      idx <- findInterval(fixed_time, fit$time)
      S_obs <- if (idx == 0 || length(fit$surv) == 0) 1 else fit$surv[idx]
      tibble(S_obs = S_obs)
    }) %>%
    ungroup()
}

#' @rdname benchmark_survival
#' @keywords internal
summarize_centers <- function(data, expected_survival, observed_survival) {
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
    mutate(O_S = S_obs * n_patients)
}

#' @rdname benchmark_survival
#' @keywords internal
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
      linewidth = 1  # MODIFIED: replaced size with linewidth
    ) +
      geom_line(
        data = lower_df,
        aes(x = precision, y = SSR, color = conf, linetype = conf),
        linewidth = 1  # MODIFIED: replaced size with linewidth
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
      geom_text(
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
  ".", "center", "status", "time", "age", "S_i", "n_patients",
  "E_S", "V_S", "O_S", "SSR", "Z", "precision", "lp", "S_obs"
))
