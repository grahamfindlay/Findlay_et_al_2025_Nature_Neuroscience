#' Statistical testing functions
#'
#' General functions for statistical hypothesis testing.

#' Get ANOVA p-value
#'
#' @param anova_result Result from anova() function
#' @return P-value from likelihood ratio test
#' @export
get_anova_pval <- function(anova_result) {
  anova_result["Pr(>Chisq)"][[1]][2]
}

#' Format post-hoc summary
#'
#' @param posthoc Post-hoc test results
#' @return Formatted summary message
#' @importFrom utils capture.output
#' @export
format_posthoc_summary <- function(posthoc) {
  # If posthoc$effect_size is a matrix, it is already in the correct format.
  # Otherwise, it is a list of lists, and we need to extract the f^2 values
  if (is.matrix(posthoc$effect_size)) {
    effect_size <- posthoc$effect_size
  } else {
    effect_size <- list()
    for (lvl in names(posthoc$effect_size)) {
      effect_size[[lvl]] <- posthoc$effect_size[[lvl]]$fsquared
    }
    effect_size <- as.matrix(effect_size)
    colnames(effect_size) <- "Cohen's f^2 analogue"
  }
  xtras <- cbind(effect_size, posthoc$ci$confint)
  msg <- capture.output(
    print(summary(posthoc$glht)),
    cat("Effect sizes and 95% family-wise confidence intervals:\n"),
    print(xtras)
  )
  msg
}

#' Test for interaction effects
#'
#' Tests for an interaction between primary and secondary covariates.
#' In the event of a significant interaction, this function performs
#' optional post-hoc tests specified by the single contrast matrix provided.
#' If you need multiple contrast matrices (e.g., one per experiment),
#' use test_interaction_with_experiment() below.
#'
#' @param dat Data frame
#' @param models List of nested models (full, reduced, null)
#' @param contrast_matrix Optional contrast matrix for post-hoc tests
#' @return List with interaction test results
#' @export
test_interaction <- function(dat, models, contrast_matrix = NULL) {
  # Test for interaction
  interaction <- list()
  interaction$anova <- stats::anova(models$full, models$reduced)
  interaction$pval <- get_anova_pval(interaction$anova)

  # If interaction is significant, estimate its effect size
  # and perform post-hoc tests
  if (interaction$pval < 0.05) {
    # Estimate interaction effect size
    interaction$effect_size <-
      subtract_ranef_get_fsquared(dat, models$full, models$reduced)

    # Perform optional post-hoc tests for nonzero slopes
    # Also get confidence intervals and effect sizes
    if (!is.null(contrast_matrix)) {
      ph <- list()
      ph$glht <- multcomp::glht(models$full, contrast_matrix)
      ph$ci <- stats::confint(ph$glht)
      ph$effect_size <- cohens_d_analogue(ph$glht, models$full)
      interaction$posthoc <- ph
    }
  }
  interaction
}

#' Test for main effects
#'
#' If the main effect is significant and a contrast matrix is provided,
#' this function performs optional post hoc tests.
#' This function does nothing to stop you from testing for a main effect
#' when such a test is incoherent (e.g., when a significant interaction
#' exists). That's on you.
#'
#' @param dat Data frame
#' @param models List of nested models (full, reduced, null)
#' @param contrast_matrix Optional contrast matrix for post-hoc tests
#' @return List with main effect test results
#' @export
test_main_effect <- function(dat, models, contrast_matrix = NULL) {
  # Test for main effect
  main_effect <- list()
  main_effect$anova <- stats::anova(models$reduced, models$null)
  main_effect$pval <- get_anova_pval(main_effect$anova)

  # If main effect is significant, estimate its effect size
  # and perform post-hoc tests
  if (main_effect$pval < 0.05) {
    # Estimate main effect size
    main_effect$effect_size <- subtract_ranef_get_fsquared(
      dat,
      models$reduced,
      models$null
    )

    # Perform optional post-hoc tests
    if (!is.null(contrast_matrix)) {
      ph <- list()
      ph$glht <- multcomp::glht(models$reduced, contrast_matrix)
      ph$ci <- stats::confint(ph$glht)
      ph$effect_size <- cohens_d_analogue(ph$glht, models$reduced)
      # Attach post-hoc tests to main effect
      main_effect$posthoc <- ph
    }
  }
  main_effect
}

#' Test interaction with experiment
#'
#' Tests for an interaction between any categorical primary covariate
#' (e.g., condition, or state) and experiment (as a secondary covariate).
#' You can also test for interactions with experiment using the more generic
#' test_interaction function above, but this function is more convenient
#' for the common case where you have a contrast matrix for each experiment.
#'
#' @param dat Data frame
#' @param models List of nested models (full, reduced, null)
#' @param novelty_matrix Contrast matrix for novelty experiment
#' @param locomotion_matrix Contrast matrix for locomotion experiment
#' @param dual_matrix Contrast matrix for dual experiment
#' @return List with interaction test results
#' @export
test_interaction_with_experiment <- function(
  dat,
  models,
  novelty_matrix,
  locomotion_matrix,
  dual_matrix
) {
  # Test for interaction
  interaction <- list()
  interaction$anova <- stats::anova(models$full, models$reduced)
  interaction$pval <- get_anova_pval(interaction$anova)

  # If interaction is significant, estimate its effect size
  # and perform post-hoc tests
  if (interaction$pval < 0.05) {
    # Estimate interaction effect size
    interaction$effect_size <- subtract_ranef_get_fsquared(
      dat,
      models$full,
      models$reduced
    )

    # Perform post-hoc tests
    ph <- list()
    # Novelty post-hoc tests, confidence intervals, and effect sizes
    ph$nod$glht <- multcomp::glht(models$full, novelty_matrix)
    ph$nod$ci <- stats::confint(ph$nod$glht)
    ph$nod$effect_size <- cohens_d_analogue(ph$nod$glht, models$full)
    # Locomotion post-hoc tests, confidence intervals, and effect sizes
    ph$cow$glht <- multcomp::glht(models$full, locomotion_matrix)
    ph$cow$ci <- stats::confint(ph$cow$glht)
    ph$cow$effect_size <- cohens_d_analogue(ph$cow$glht, models$full)
    # Dual post-hoc tests, confidence intervals, and effect sizes
    ph$ctn$glht <- multcomp::glht(models$full, dual_matrix)
    ph$ctn$ci <- stats::confint(ph$ctn$glht)
    ph$ctn$effect_size <- cohens_d_analogue(ph$ctn$glht, models$full)
    # Attach post-hoc tests to interaction
    interaction$posthoc <- ph
  }
  interaction
}
