#' Benchmark Survival Analysis
#'
#' Performs survival benchmarking across centers using standardized survival ratios and funnel plots.
#'
#' @param data A data frame containing survival data.
#' @param center_col Name of the column with center identifiers (character).
#' @param status_col Name of the column with event status (0 = censored, 1 = event).
#' @param time_col Name of the column with survival times (numeric).
#' @param fixed_time Time point for survival analysis (e.g., 36 for 3 years).
#' @param numeric_covariates Optional vector of column names for numeric covariates used for risk adjustment.
#' @param factor_covariates Optional vector of column names for categorical covariates used for risk adjustment.
#' @param time_dep_covariates Optional vector of numeric covariate names that violate the PH assumption. These must be included in \code{numeric_covariates}. For each, an interaction with \code{log(time)} is included in the model.
#' @param strata_covariates Optional vector of column names for stratification via \code{strata()}.
#' @param highlight_center Optional center identifier to highlight in the plot.
#' @param conf_levels Numeric vector of confidence levels for funnel plot limits (e.g., \code{c(0.8, 0.95)}).
#' @param include_legend Logical; include legend in the plot (default \code{TRUE}).
#' @return A \code{ggplot} object representing the funnel plot.
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
  # Input validation
  if (!is.data.frame(data)) stop("`data` must be a data frame")
  if (!is.character(center_col) || !is.character(status_col) || !is.character(time_col)) {
    stop("`center_col`, `status_col`, and `time_col` must be character strings")
  }
  if (!is.numeric(fixed_time) || fixed_time <= 0) stop("`fixed_time` must be a positive numeric value")
  if (!all(conf_levels > 0 & conf_levels < 1)) stop("`conf_levels` must be between 0 and 1")

  # Prepare data
  prepared_data <- prepare_data(data, center_col, status_col, time_col,
                                covariates = c(numeric_covariates, factor_covariates),
                                factor_covariates = factor_covariates,
                                strata_covariates = strata_covariates)

  # Fit Cox model
  cox_model <- fit_cox_model(prepared_data, c(numeric_covariates, factor_covariates),
                             strata_covariates, time_dep_covariates, fixed_time)

  # Check PH assumption
  ph_test <- survival::cox.zph(cox_model)
  message("\nProportional Hazards (PH) Assumption Test:")
  print(ph_test)
  if (all(ph_test$table[, "p"] > 0.05)) {
    message("✅ The Proportional Hazards (PH) assumption is fulfilled.")
  } else {
    warning("⚠️ PH assumption may be violated. Consider modeling violating covariates as time-dependent.")
    message("Covariates violating PH assumption: ",
            paste(rownames(ph_test$table)[ph_test$table[, "p"] <= 0.05], collapse = ", "))
  }

  # Compute survival metrics and plot
  expected_survival <- compute_expected_survival(prepared_data, cox_model, fixed_time)
  observed_survival <- compute_observed_survival(prepared_data, fixed_time)
  center_summary <- summarize_centers(prepared_data, expected_survival, observed_survival)
  metrics <- compute_metrics(center_summary)
  funnel_limits <- generate_funnel_limits(metrics, conf_levels)
  plot <- create_funnel_plot(metrics, funnel_limits, highlight_center, conf_levels, include_legend)

  # Model fit summary
  cox_fit_summary <- summary(cox_model)
  message("\nModel Fit Assessment:")
  message("✔️ Concordance Index (C-statistic): ", round(cox_fit_summary$concordance[1], 3))
  message("✔️ Likelihood Ratio Test p-value: ", format.pval(cox_fit_summary$logtest["pvalue"], digits = 3))
  message("✔️ Wald Test p-value: ", format.pval(cox_fit_summary$waldtest["pvalue"], digits = 3))
  message("✔️ Score (Log-Rank) Test p-value: ", format.pval(cox_fit_summary$sctest["pvalue"], digits = 3))

  return(plot)
}

#' @rdname benchmark_survival
#' @keywords internal
prepare_data <- function(data, center_col, status_col, time_col, covariates,
                         factor_covariates, strata_covariates = NULL) {
  required_cols <- unique(c(center_col, status_col, time_col, covariates, strata_covariates))
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

  data <- data %>%
    rename(center = !!sym(center_col),
           status = !!sym(status_col),
           time = !!sym(time_col)) %>%
    mutate(status = as.integer(status),
           time = as.numeric(time)) %>%
    filter(!is.na(time), !is.na(status), time > 0)

  if (!is.null(factor_covariates)) {
    data <- data %>% mutate(across(all_of(factor_covariates), as.factor))
  }
  if (!is.null(strata_covariates)) {
    data <- data %>% mutate(across(all_of(strata_covariates), as.factor))
  }

  return(data)
}

#' @rdname benchmark_survival
#' @keywords internal
fit_cox_model <- function(data, covariates, strata_covariates = NULL,
                          time_dep_covariates = NULL, fixed_time) {
  if (!is.null(time_dep_covariates)) {
    if (!all(time_dep_covariates %in% covariates)) {
      stop("All `time_dep_covariates` must be in `covariates`")
    }
  }

  # Build formula
  risk_terms <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else "1"
  if (!is.null(time_dep_covariates)) {
    td_terms <- paste0("tt(", time_dep_covariates, ")", collapse = " + ")
    risk_terms <- paste(risk_terms, td_terms, sep = " + ")
  }
  if (!is.null(strata_covariates)) {
    strata_terms <- paste0("strata(", strata_covariates, ")", collapse = " + ")
    formula_str <- paste("Surv(time, status) ~", risk_terms, "+", strata_terms)
  } else {
    formula_str <- paste("Surv(time, status) ~", risk_terms)
  }

  # Fit model
  tt <- function(x, t, ...) x * log(t)
  survival::coxph(as.formula(formula_str), data = data, tt = tt)
}

#' @rdname benchmark_survival
#' @keywords internal
compute_expected_survival <- function(data, cox_model, fixed_time) {
  bh <- survival::basehaz(cox_model, centered = FALSE)
  max_time <- max(bh$time, na.rm = TRUE)
  if (max_time < fixed_time) {
    warning("Maximum follow-up time (", round(max_time, 1),
            ") is less than `fixed_time` (", fixed_time, "). Results may be unreliable.")
  }
  H_fixed <- max(bh$hazard[bh$time <= fixed_time], 0, na.rm = TRUE)
  lp <- predict(cox_model, newdata = data, type = "lp")

  data %>% mutate(S_i = exp(-H_fixed * exp(lp)))
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
    summarise(E_S = sum(S_i, na.rm = TRUE),
              V_S = sum(S_i * (1 - S_i), na.rm = TRUE),
              .groups = "drop")

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
    mutate(SSR = O_S / E_S,
           Z = (O_S - E_S) / sqrt(V_S),
           precision = E_S^2 / V_S) %>%
    filter(is.finite(SSR), is.finite(precision), precision > 0)
}

#' @rdname benchmark_survival
#' @keywords internal
generate_funnel_limits <- function(metrics, conf_levels) {
  precision_range <- range(metrics$precision, na.rm = TRUE)
  precision_seq <- seq(max(precision_range[1], 0.01), precision_range[2], length.out = 100)
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
  color_map <- setNames(c("blue", "red"), paste0(conf_levels * 100, "%"))
  linetype_map <- setNames(c("dashed", "dotted"), paste0(conf_levels * 100, "%"))
  y_min <- min(metrics$SSR, na.rm = TRUE) - 0.1

  gg <- ggplot() +
    geom_point(data = metrics,
               aes(x = precision, y = SSR, fill = center),
               shape = 21, color = "black", size = 3, alpha = 0.7) +
    geom_hline(yintercept = 1, color = "grey50", linetype = 2) +
    theme_minimal(base_size = 14) +
    labs(title = "Survival Outcomes Benchmarking",
         x = "Precision (E^2/V)",
         y = "Standardized Survival Ratio (O/E)",
         fill = "Center",
         color = "Confidence",
         linetype = "Confidence") +
    scale_fill_discrete() +
    scale_color_manual(values = color_map, drop = FALSE) +
    scale_linetype_manual(values = linetype_map, drop = FALSE) +
    coord_cartesian(ylim = c(y_min, NA)) +
    theme(legend.position = if (include_legend) "right" else "none")

  for (conf in conf_levels) {
    conf_label <- paste0(conf * 100, "%")
    upper_df <- tibble(precision = limits$precision,
                       SSR = limits[[paste0("upper_", conf)]],
                       conf = conf_label)
    lower_df <- tibble(precision = limits$precision,
                       SSR = limits[[paste0("lower_", conf)]],
                       conf = conf_label)

    gg <- gg +
      geom_line(data = upper_df, aes(x = precision, y = SSR, color = conf, linetype = conf), size = 1) +
      geom_line(data = lower_df, aes(x = precision, y = SSR, color = conf, linetype = conf), size = 1)
  }

  if (!is.null(highlight)) {
    highlighted_data <- filter(metrics, center == highlight)
    if (nrow(highlighted_data) > 0) {
      gg <- gg +
        geom_point(data = highlighted_data,
                   aes(x = precision, y = SSR),
                   shape = 21, fill = "red", color = "black", size = 5, stroke = 1.3) +
        geom_text(data = highlighted_data,
                  aes(x = precision, y = SSR, label = center),
                  vjust = -1, color = "red", size = 5, fontface = "bold")
    }
  }

  return(gg)
}

utils::globalVariables(c("center", "status", "time", "S_i", "n_patients",
                         "E_S", "V_S", "O_S", "SSR", "Z", "precision", "lp", "S_obs"))
