#' Homeostasis testing functions
#'
#' Functions specific to testing homeostasis hypotheses.

#' Test homeostasis hypothesis
#'
#' Tests for homeostasis effects across sleep/wake conditions and experiments.
#'
#' @param response_variable Name of the response variable to test
#' @return List with homeostasis test results
#' @export
test_homeostasis <- function(response_variable) {
  r <- list()
  r$data <- load_measures_for_homeostasis_testing()
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
    homeostasis_novelty_matrix,
    homeostasis_locomotion_matrix,
    homeostasis_dual_matrix
  )
  if (r$interaction$pval >= 0.05) {
    r$main_effect <- test_main_effect(
      r$data,
      r$models,
      homeostasis_main_effect_matrix
    )
  }
  r
}
