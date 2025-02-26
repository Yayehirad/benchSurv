#' Benchmark Survival Analysis
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
#' @importFrom survival coxph basehaz survfit Surv
#' @importFrom ggplot2 ggplot geom_point geom_line labs scale_color_manual theme_minimal coord_cartesian geom_hline scale_linetype_manual aes theme scale_fill_discrete
#' @importFrom colorspace rainbow_hcl
#' @importFrom stats coef qnorm setNames
#' @importFrom utils tail globalVariables
#' @importFrom rlang sym .data
#' @importFrom tibble tibble
#' @examples
#' set.seed(123)
#' data <- data.frame(
#'   SiteShortName = rep(c("SiteA", "SiteB", "SiteC", "SiteD", "SiteE"), each = 100),
#'   Death = sample(c(0, 1), 100, replace = TRUE),
#'   DTime = sample(10:100, 100, replace = TRUE),
#'   Age = sample(30:80, 100, replace = TRUE)
#' )
#' plot <- benchmark_survival(data, "SiteShortName", "Death", "DTime", "Age", 36, "SiteA")
benchmark_survival <- function(data, center_col, status_col, time_col, age_col, fixed_time,
                               highlight_center = NULL, conf_levels = c(0.8, 0.95)) {
  if (!is.data.frame(data)) stop("data must be a data frame") # Data validation
  if (!is.numeric(fixed_time) || fixed_time <= 0) stop("fixed_time must be positive numeric")

  # Processing pipeline
  prepared_data <- prepare_data(data, center_col, status_col, time_col, age_col)
  cox_model <- fit_cox_model(prepared_data)
  expected_survival <- compute_expected_survival(prepared_data, cox_model, fixed_time)
  observed_survival <- compute_observed_survival(prepared_data, fixed_time)
  center_summary <- summarize_centers(prepared_data, expected_survival, observed_survival)
  metrics <- compute_metrics(center_summary)
  funnel_limits <- generate_funnel_limits(metrics, conf_levels)
  create_funnel_plot(metrics, funnel_limits, highlight_center, conf_levels)
}

# Internal utility functions --------------------------------------------------

#' @rdname benchmark_survival
#' @keywords internal
prepare_data <- function(data, center_col, status_col, time_col, age_col) {
  required_cols <- c(center_col, status_col, time_col, age_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

  data %>%
    rename(
      center = !!sym(center_col),
      status = !!sym(status_col),
      time   = !!sym(time_col),
      age    = !!sym(age_col)
    ) %>%
    mutate(
      age = as.numeric(age),
      status = as.integer(status),
      time = as.numeric(time)
    ) %>%
    filter(
      !is.na(age),
      !is.na(time),
      !is.na(status),
      time > 0
    )
}

#' @rdname benchmark_survival
#' @keywords internal
fit_cox_model <- function(data) {
  survival::coxph(Surv(time, status) ~ age, data = data)
}

#' @rdname benchmark_survival
#' @keywords internal
compute_expected_survival <- function(data, cox_model, fixed_time) {
  bh <- survival::basehaz(cox_model, centered = FALSE)
  H_fixed <- max(bh$hazard[bh$time <= fixed_time], 0)

  data %>%
    mutate(
      lp = coef(cox_model) * age,
      S_i = exp(-H_fixed * exp(lp))
    )
}

#' @rdname benchmark_survival
#' @keywords internal
compute_observed_survival <- function(data, fixed_time) {
  data %>%
    group_by(center) %>%
    group_modify(~ {
      # Extend the survival curve so it covers fixed_time even if there are no events after the last observation
      fit <- survival::survfit(Surv(time, status) ~ 1, data = .x)
      idx <- findInterval(fixed_time, fit$time)
      S_obs <- if (idx == 0) 1 else fit$surv[idx]
      tibble::tibble(S_obs = S_obs)
    }) %>%
    ungroup()
}

#' @rdname benchmark_survival
#' @keywords internal
summarize_centers <- function(data, expected_survival, observed_survival) {
  # Aggregate patient-level expected survival and then join with observed survival
  patient_agg <- expected_survival %>%
    group_by(center) %>%
    summarise(
      E_S = sum(S_i),
      V_S = sum(S_i * (1 - S_i)),
      .groups = "drop"
    )

  data %>%
    group_by(center) %>%
    summarise(
      n_patients = dplyr::n(),
      .groups = "drop"
    ) %>%
    left_join(patient_agg, by = "center") %>%
    left_join(observed_survival, by = "center") %>%
    mutate(
      O_S = S_obs * n_patients
    )
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
create_funnel_plot <- function(metrics, limits, highlight = NULL, conf_levels) {
  # Define color and linetype mappings for confidence levels
  color_map <- c("80%" = "blue", "95%" = "red")
  linetype_map <- c("80%" = "dashed", "95%" = "dotted")

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
    # Horizontal reference line
    geom_hline(yintercept = 1, color = "grey50", linetype = 2) +
    theme_minimal(base_size = 14) +
    labs(
      title  = "Survival Outcomes Benchmarking",
      x      = "Precision (E^2/V)",
      y      = "Standardized Survival Ratio",
      fill   = "Center",       # legend title for center fills
      color  = "Confidence",   # legend title for confidence lines
      linetype = "Confidence"
    ) +
    # Apply discrete color/fill scales for centers
    ggplot2::scale_fill_discrete() +
    ggplot2::scale_color_manual(values = color_map, drop = FALSE) +
    scale_linetype_manual(values = linetype_map, drop = FALSE) +
    # Y-axis from 0 to max
    coord_cartesian(ylim = c(0, NA)) +
    # Place legend in top-right corner
    theme(
      legend.position = c(1, 1),        # top-right corner
      legend.justification = c(1, 1)    # anchor the legend at that corner
    )

  # Plot the upper and lower lines for each confidence level in separate geoms
  for (conf in conf_levels) {
    conf_label <- paste0(conf * 100, "%")

    # Upper line
    upper_df <- tibble::tibble(
      precision = limits$precision,
      SSR       = limits[[paste0("upper_", conf)]],
      conf      = conf_label
    )

    # Lower line
    lower_df <- tibble::tibble(
      precision = limits$precision,
      SSR       = limits[[paste0("lower_", conf)]],
      conf      = conf_label
    )

    gg <- gg + geom_line(
      data = upper_df,
      aes(x = precision, y = SSR, color = conf, linetype = conf),
      size = 1
    )
    gg <- gg + geom_line(
      data = lower_df,
      aes(x = precision, y = SSR, color = conf, linetype = conf),
      size = 1
    )
  }

  # Highlight a center if specified
  if (!is.null(highlight)) {
    gg <- gg +
      geom_point(
        data = dplyr::filter(metrics, center == highlight),
        aes(x = precision, y = SSR),
        shape = 21,
        fill = "yellow",
        color = "black",
        size = 5,
        stroke = 1.3,
        show.legend = FALSE
      )
  }

  gg
}


# Global variable declarations to avoid R CMD check notes
utils::globalVariables(c(
  ".", "center", "status", "time", "age", "S_i", "n_patients",
  "E_S", "V_S", "O_S", "SSR", "Z", "precision", "lp", "S_obs"
))
