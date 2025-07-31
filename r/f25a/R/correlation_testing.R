#' Correlation testing functions
#'
#' Functions for correlation analyses between conditions and variables.

#' Get correlation data
#'
#' Transform data for correlation analysis between two condition/variable pairs.
#' If x_condition = "EXT.Wake", x_variable = "Theta",
#' y_condition = "REC.NREM", y_variable = "SPWs"...
#'
#' Transform data d of the form:
#'
#' subject experiment condition THETA SPWs
#' 1       A          EXT.Wake   10    20
#' 1       A          REC.NREM   15    25
#' 1       B          EXT.Wake   20    30
#' 1       B          REC.NREM   25    35
#'
#' into dxy of the form:
#'
#' subject experiment `EXT.Wake.Theta` `REC.NREM.SPWs`
#' 1       A          10              25
#' 1       B          20              35
#' 2       A          ...             ...
#'
#' where xvar = "`EXT.Wake.Theta`" and yvar = "`REC.NREM.SPWs`"
#'
#' @param d Data frame in long format
#' @param x_condition Condition name for x variable
#' @param x_variable Variable name for x
#' @param y_condition Condition name for y variable
#' @param y_variable Variable name for y
#' @return List with merged data frame and variable names
#' @export
get_correlation_data <- function(
  d,
  x_condition,
  x_variable,
  y_condition,
  y_variable
) {
  xvar <- paste(x_condition, x_variable, sep = ".")
  dx <- d[d$condition == x_condition, ]
  dx[xvar] <- dx[x_variable]
  dx <- dx[, c("subject", "experiment", xvar)]

  yvar <- paste(y_condition, y_variable, sep = ".")
  dy <- d[d$condition == y_condition, ]
  dy[yvar] <- dy[y_variable]
  dy <- dy[, c("subject", "experiment", yvar)]

  dxy <- merge(dx, dy, by = c("subject", "experiment"))
  list(data = dxy, xvar = xvar, yvar = yvar)
}

#' Fit dummy models for f-squared calculation
#'
#' If you have covariate X (hippocampal theta-delta ratio),
#' categorical variable Z (experiment) with levels A (NOD), B (COW), C (CTN).
#'
#' The full model y ~ X * Z will have terms
#' \code{[intercept, zB, zC, X, zB:X, zC:X]}.
#' To test for an effect for A (Novelty), you test if \code{X = 0}.
#' To test for an effect for B (Locomotion), you test if \code{X + zB:X = 0}.
#' To test for an effect for C (Locomotion + Novelty), you test if
#' \code{X + zC:X = 0}.
#'
#' @param d Data frame
#' @param response_var Response variable name
#' @param primary_covariate Primary covariate name
#' @param secondary_covariate Secondary covariate name
#' @param random_factor Random factor name
#' @return List with dummy-coded data and models
#' @export
fit_dummy_models_for_fsquared <- function(
  d,
  response_var,
  primary_covariate,
  secondary_covariate,
  random_factor
) {
  dummy_d <- fastDummies::dummy_cols(
    d,
    select_columns = secondary_covariate,
    omit_colname_prefix = TRUE
  )
  dummy_vars <- fastDummies::dummy_cols(
    d,
    select_columns = secondary_covariate,
    omit_colname_prefix = TRUE,
    return_generated_variables = TRUE
  )

  random_intercept <- paste0("(1 | ", random_factor, ")")

  # Create formula string for full dummy model
  # The "0" forces the intercept to be zero, which allows the dummy variables
  # to appear in the model, which would otherwise be overdetermined.
  simple_fixefs <- paste(dummy_vars, collapse = " + ")
  intractn_terms <- paste0(dummy_vars, ":", primary_covariate)
  intractn_fixefs <- paste(intractn_terms, collapse = " + ")
  dummy_formula <- stats::reformulate(
    c(
      paste("0", "+", simple_fixefs, "+", intractn_fixefs),
      random_intercept
    ),
    response = response_var
  )

  # Fit full dummy model
  dummy_model <- lme4::lmer(dummy_formula, data = dummy_d, REML = FALSE)

  # Create a list to store models with each interaction term removed
  reduced_models <- list()

  # For each interaction term, create a model without that term
  for (i in seq_along(intractn_terms)) {
    term <- intractn_terms[i]
    dummy <- dummy_vars[i]

    # Remove the current interaction term from the formula
    remaining_terms <- intractn_terms[intractn_terms != term]
    remaining_fixefs <- paste(remaining_terms, collapse = " + ")

    # Create new formula without this interaction term
    reduced_formula <- stats::reformulate(
      c(
        paste("0", "+", simple_fixefs, "+", remaining_fixefs),
        random_intercept
      ),
      response = response_var
    )

    # Fit reduced model and store it
    reduced_models[[dummy]] <- lme4::lmer(
      reduced_formula,
      data = dummy_d,
      REML = FALSE
    )
  }

  list(
    data = dummy_d,
    models = list(
      full = dummy_model,
      reduced = reduced_models
    )
  )
}

#' Get f-squared for interaction post-hocs
#'
#' Get Cohen's f-squared for posthoc tests, when Cohen's d is not appropriate.
#' The issue is that Cohen's d estimates may depend on the scale of the inputs.
#' You can change the scale of the inputs, making the estimates larger or
#' smaller arbitrarily, and the residual error won't change (linear
#' transformation of a linear model). So that notion of effect size doesn't
#' work, but something like Cohen's f-squared does.
#'
#' @param d Data frame
#' @param response_var Response variable name
#' @param primary_covariate Primary covariate name
#' @param secondary_covariate Secondary covariate name
#' @param random_factor Random factor name
#' @return List of effect sizes by level
#' @export
get_fsquared_for_interaction_posthocs <- function(
  d,
  response_var,
  primary_covariate,
  secondary_covariate,
  random_factor
) {
  dummy <- fit_dummy_models_for_fsquared(
    d,
    response_var,
    primary_covariate,
    secondary_covariate,
    random_factor
  )

  effect_sizes <- list()
  for (lvl in names(dummy$models$reduced)) {
    effect_sizes[[lvl]] <- subtract_ranef_get_fsquared(
      dummy$data,
      dummy$models$full,
      dummy$models$reduced[[lvl]]
    )
  }
  effect_sizes
}

#' Test condition correlation
#'
#' Test correlation between variables from different conditions.
#'
#' @param x_condition Condition name for x variable
#' @param x_variable Variable name for x
#' @param y_condition Condition name for y variable
#' @param y_variable Variable name for y
#' @return List with correlation test results
#' @export
test_condition_correlation <- function(
  x_condition,
  x_variable,
  y_condition,
  y_variable
) {
  measures <- load_measures()
  dxy <- get_correlation_data(
    measures,
    x_condition,
    x_variable,
    y_condition,
    y_variable
  )

  r <- list() # Hold results
  r$data <- dxy$data
  r$xvar <- dxy$xvar
  r$yvar <- dxy$yvar
  r$models <- fit_nested_models(
    dxy$data,
    dxy$yvar,
    dxy$xvar,
    "experiment",
    "subject"
  )

  r$interaction <- test_interaction(
    r$data,
    r$models,
    correlation_interaction_matrix
  )
  r$interaction$posthoc$effect_size <- get_fsquared_for_interaction_posthocs(
    r$data,
    dxy$yvar,
    dxy$xvar,
    "experiment",
    "subject"
  ) # Overwrite faulty Cohen's d analogues with f-squared
  assertthat::assert_that(
    all(
      rownames(correlation_interaction_matrix) ==
        names(r$interaction$posthoc$effect_size)
    ),
    msg = "Contrast matrix rows and posthoc levels do not exactly match!"
  ) # Not strictly necessary, but nice for later formatting of outputs.
  if (r$interaction$pval >= 0.05) {
    r$main_effect <- test_main_effect(r$data, r$models)
  }
  r
}

#' Test condition difference correlation
#'
#' Test correlation between variables using condition difference data.
#'
#' @param x_variable Variable name for x
#' @param y_variable Variable name for y
#' @return List with correlation test results
#' @export
test_condition_diff_corr <- function(x_variable, y_variable) {
  r <- list()
  r$data <- load_condition_differences()
  r$models <- fit_nested_models(
    r$data,
    y_variable,
    x_variable,
    "experiment",
    "subject"
  )

  r$interaction <- test_interaction(
    r$data,
    r$models,
    correlation_interaction_matrix
  )
  r$interaction$posthoc$effect_size <- get_fsquared_for_interaction_posthocs(
    r$data,
    y_variable,
    x_variable,
    "experiment",
    "subject"
  ) # Overwrite faulty Cohen's d analogues with f-squared
  assertthat::assert_that(
    all(
      rownames(correlation_interaction_matrix) ==
        names(r$interaction$posthoc$effect_size)
    ),
    msg = "Contrast matrix rows and posthoc levels do not exactly match!"
  ) # Not strictly necessary, but nice for later formatting of outputs.
  if (r$interaction$pval >= 0.05) {
    r$main_effect <- test_main_effect(r$data, r$models)
  }
  r
}
