#' Other specialized testing functions
#'
#' Functions for extended wake change testing, occupancy testing, and other specialized analyses.

#' Test extended wake change
#'
#' Tests for changes during extended wake periods across experiments.
#'
#' @param response_variable Name of the response variable to test
#' @return List with extended wake change test results
#' @export
test_ewk_change <- function(response_variable) {
  r <- list()
  r$data <- load_measures_for_ewk_change_testing()
  r$models <- fit_nested_models(
    r$data,
    response_variable,
    "condition",
    "experiment",
    "subject"
  )

  r$interaction <- test_interaction_with_experiment(
    r$data,
    r$models,
    ewk_change_novelty_matrix,
    ewk_change_locomotion_matrix,
    ewk_change_dual_matrix
  )
  if (r$interaction$pval >= 0.05) {
    r$main_effect <- test_main_effect(
      r$data,
      r$models,
      ewk_change_main_effect_matrix
    )
  }
  r
}

#' Test single day occupancy
#'
#' Tests occupancy differences between experiments for a single day.
#'
#' @param day Day to test ("Baseline" or "Recovery")
#' @return List with occupancy test results
#' @export
test_single_day_occupancy <- function(day) {
  d <- load_occupancy_data()

  r <- list()
  r$data <- droplevels(d[d$day == day, ])
  r$models <- fit_nested_models(
    r$data,
    "fractional_occupancy",
    "experiment",
    "state",
    "subject"
  )

  r$interaction <- test_interaction(
    r$data,
    r$models,
    day_occupancy_contrast_matrix
  )
  if (r$interaction$pval >= 0.05) {
    r$main_effect <- test_main_effect(r$data, r$models)
    assertthat::assert_that(
      r$main_effect$pval >= 0.05,
      msg = "Contrast matrix for main effect posthocs not implemented!"
    )
  }
  r
}

#' Test single experiment occupancy
#'
#' Tests occupancy differences between days for a single experiment.
#'
#' @param experiment Experiment to test ("Novelty", "Locomotion", or "Dual")
#' @return List with occupancy test results
#' @export
test_single_experiment_occupancy <- function(experiment) {
  d <- load_occupancy_data()

  r <- list()
  r$data <- droplevels(d[d$experiment == experiment, ])
  r$models <- fit_nested_models(
    r$data,
    "fractional_occupancy",
    "day",
    "state",
    "subject"
  )

  r$interaction <- test_interaction(
    r$data,
    r$models,
    exp_occupancy_contrast_matrix
  )
  if (r$interaction$pval >= 0.05) {
    r$main_effect <- test_main_effect(r$data, r$models)
    assertthat::assert_that(
      r$main_effect$pval >= 0.05,
      msg = "Contrast matrix for main effect posthocs not implemented!"
    )
  }
  r
}
