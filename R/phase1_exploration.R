# R/phase1_exploration.R
# Phase 1: goal declaration, data inspection, and off-ramp assessment.
# Called interactively from templates/phase1_exploration.qmd.
# DESIGN.md §§4, 5.

run_phase1 <- function(wf, data = NULL, outcome_var = NULL, outcome_type = NULL,
                       goal = NULL, simulated = FALSE) {
  stopifnot(inherits(wf, "wf_state"))

  # --- Simulated-data fast path ---
  if (simulated) {
    if (is.null(goal)) goal <- "model comparison on simulated data with known ground truth"
    wf$declared_goal  <- goal
    wf$n_observations <- NA
    cat("=== Phase 1: Simulated Data ===\n")
    cat("Goal:", goal, "\n")
    cat("Off-ramp assessment skipped (ground truth is known).\n\n")
    wf$audit_trail <- append(wf$audit_trail, list(list(
      phase     = 1,
      action    = "bayesian_selected",
      timestamp = Sys.time(),
      method    = "full_stan",
      goal      = goal
    )))
    export_context(wf)
    return(invisible(wf))
  }

  stopifnot(is.data.frame(data))
  stopifnot(is.character(outcome_var), length(outcome_var) == 1)
  stopifnot(outcome_type %in% c("binary", "count", "continuous"))
  stopifnot(is.character(goal), length(goal) == 1)

  # --- Store metadata on wf ---
  wf$declared_goal  <- goal
  wf$n_observations <- nrow(data)
  if (outcome_type == "binary" && outcome_var %in% names(data)) {
    y_raw      <- data[[outcome_var]]
    wf$event_rate <- mean(as.integer(y_raw), na.rm = TRUE)
  }

  # --- Exploratory plot ---
  cat("=== Phase 1: Data Inspection ===\n")
  cat("Outcome variable:", outcome_var, "\n")
  cat("Outcome type:    ", outcome_type, "\n")
  cat("n observations:  ", wf$n_observations, "\n")
  if (!is.null(wf$event_rate) && !is.na(wf$event_rate))
    cat("Event rate:      ", sprintf("%.1f%%\n", wf$event_rate * 100))
  cat("Declared goal:   ", goal, "\n\n")

  p <- if (outcome_type %in% c("count", "continuous")) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[outcome_var]])) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
      ggplot2::labs(title = paste("Distribution of", outcome_var),
                    x = outcome_var, y = "Count") +
      ggplot2::theme_minimal()
  } else {
    tbl <- table(data[[outcome_var]])
    bar_df <- data.frame(value = names(tbl), count = as.integer(tbl))
    ggplot2::ggplot(bar_df, ggplot2::aes(x = .data[["value"]],
                                          y = .data[["count"]])) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::labs(title = paste("Distribution of", outcome_var),
                    x = outcome_var, y = "Count") +
      ggplot2::theme_minimal()
  }
  print(p)

  # --- Off-ramp assessment ---
  offramp_result <- assess_offramps(
    data         = data,
    outcome_var  = outcome_var,
    outcome_type = outcome_type,
    goal         = goal,
    n            = wf$n_observations
  )
  print_offramps(offramp_result)

  # --- User choice ---
  choices     <- seq_along(offramp_result$alternatives)
  choice_strs <- paste(choices, collapse = "/")
  ack <- readline(sprintf(
    "Select an analysis path [%s], or press Enter for Full Bayesian: ",
    choice_strs
  ))

  selected_idx <- suppressWarnings(as.integer(trimws(ack)))
  if (is.na(selected_idx) || selected_idx < 1 ||
      selected_idx > length(offramp_result$alternatives)) {
    # Default: Full Bayesian
    action <- "bayesian_selected"
    method <- "full_stan"
    cat("Full Bayesian path selected.\n")
  } else {
    method <- offramp_result$alternatives[[selected_idx]]$method
    if (method == "full_stan") {
      action <- "bayesian_selected"
      cat("Full Bayesian path selected.\n")
    } else {
      action <- "offramp_choice"
      cat("Off-ramp selected:", offramp_result$alternatives[[selected_idx]]$label, "\n")
    }
  }

  # --- Log to audit trail ---
  wf$audit_trail <- append(wf$audit_trail, list(list(
    phase     = 1,
    action    = action,
    timestamp = Sys.time(),
    method    = method,
    goal      = goal
  )))

  export_context(wf)
  invisible(wf)
}
