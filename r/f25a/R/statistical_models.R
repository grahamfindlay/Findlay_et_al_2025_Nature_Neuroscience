#' Statistical modeling functions
#'
#' Core functions for fitting nested models and effect size calculations.

#' Fit nested statistical models
#'
#' Fits full, reduced, and null models for nested model comparisons.
#'
#' @param d Data frame
#' @param response_var Response variable name
#' @param primary_covariate Primary covariate name
#' @param secondary_covariate Secondary covariate name
#' @param random_factor Random factor name
#' @return List of fitted models (full, reduced, null)
#' @export
fit_nested_models <- function(
  d,
  response_var,
  primary_covariate,
  secondary_covariate,
  random_factor
) {
  random_intercept <- paste0("(1 | ", random_factor, ")")

  # Create formula strings
  f1 <- stats::reformulate(
    c(paste(primary_covariate, "*", secondary_covariate), random_intercept),
    response = response_var
  ) # Full model

  f2 <- stats::reformulate(
    c(primary_covariate, secondary_covariate, random_intercept),
    response = response_var
  ) # Reduced model

  f3 <- stats::reformulate(
    c(secondary_covariate, random_intercept),
    response = response_var
  ) # Null model

  # Fit models
  m1 <- lme4::lmer(f1, data = d, REML = FALSE) # Full model
  m2 <- lme4::lmer(f2, data = d, REML = FALSE) # Reduced model
  m3 <- lme4::lmer(f3, data = d, REML = FALSE) # Null model

  list(full = m1, reduced = m2, null = m3)
}

#' Fit basic statistical models
#'
#' Fits basic full and null models for simpler comparisons.
#'
#' @param d Data frame
#' @param response_var Response variable name
#' @param primary_covariate Primary covariate name
#' @param random_factor Random factor name
#' @return List of fitted models (full, null)
#' @export
fit_basic_models <- function(
  d,
  response_var,
  primary_covariate,
  random_factor
) {
  primary_covariate <- paste0("`", primary_covariate, "`")
  response_var <- paste0("`", response_var, "`")
  random_factor <- paste0("`", random_factor, "`")
  random_intercept <- paste0("(1 | ", random_factor, ")")

  # Create formula strings
  f1 <- stats::reformulate(
    c(primary_covariate, random_intercept),
    response = response_var
  ) # Full model

  f2 <- stats::reformulate(
    c(random_intercept),
    response = response_var
  ) # Null model

  # Fit models
  m1 <- lme4::lmer(f1, data = d, REML = FALSE) # Full model
  m2 <- lme4::lmer(f2, data = d, REML = FALSE) # Null model

  list(full = m1, null = m2)
}

#' Subtract simple (non-crossed, non-nested, single-level) random effects
#'
#' @param d The data with random effects included, used to fit model m
#' @param m The model with random effects structure that must be removed
#' @return Data frame with random effects subtracted
#' @importFrom utils tail
#' @export
subtract_simple_random_effects <- function(
  d, # The data with random effects included, used to fit model m.
  m # The model with random effects structure that must be removed.
) {
  fixed <- d # Make a copy of the data

  # Extract model components. If the model is:
  #   Cortical SWA ~ Condition * Experiment + (1 | Subject)
  # then Cortical SWA is the response variable, Condition and Experiment are
  # the primary and secondary covariates, and Subject is the random factor.
  random_effects <- lme4::ranef(m)
  random_factor <- names(random_effects)
  response_var <- all.vars(stats::formula(m))[1]
  assertthat::assert_that(
    random_factor == tail(all.vars(stats::formula(m)), n = 1)
  )

  # For each random factor level (e.g., for each subject),
  # subtract the random effect from the response variable.
  random_factor_levels <- levels(fixed[[random_factor]])
  for (fctr in random_factor_levels) {
    intercept <- random_effects[[random_factor]][fctr, 1]
    is_fctr <- (fixed[[random_factor]] == fctr)
    fixed[is_fctr, response_var] <- fixed[is_fctr, response_var] - intercept
  }
  fixed
}

#' Subtract nested random effects from a data frame
#'
#' @param d The data with random effects included, used to fit model m
#' @param m The model with random effects structure that must be removed
#' @return Data frame with random effects subtracted
#' @export
subtract_nested_random_effects <- function(d, m) {
  fixed <- d # Make a copy of the data

  random_effects <- lme4::ranef(m)
  re_names <- names(random_effects)
  nested_name <- re_names[grepl(":", re_names)]
  if (length(nested_name) == 0) {
    stop("No nested random effects found in model")
  }
  if (length(nested_name) > 1) {
    stop(
      "Multiple nested random effects found. ",
      "This function supports only one nested structure."
    )
  }

  # Parse the nested name to get variable names (format: "inner:outer")
  nested_parts <- strsplit(nested_name, ":")[[1]]
  inner_var <- nested_parts[1]
  outer_var <- nested_parts[2]

  # Find the outer-level random effect name
  outer_name <- re_names[re_names == outer_var]
  if (length(outer_name) == 0) {
    stop(
      "Outer grouping variable ",
      outer_var,
      " not found in random effects"
    )
  }

  # Subtract outer-level random effects
  response_var <- all.vars(stats::formula(m))[1]
  outer_effects <- random_effects[[outer_name]]
  for (outer_id in rownames(outer_effects)) {
    outer_intercept <- outer_effects[outer_id, 1]
    is_outer <- (fixed[[outer_var]] == outer_id)
    fixed[is_outer, response_var] <-
      fixed[is_outer, response_var] - outer_intercept
  }

  # Subtract inner-level random effects (nested within outer)
  nested_effects <- random_effects[[nested_name]]
  for (nested_id in rownames(nested_effects)) {
    nested_intercept <- nested_effects[nested_id, 1]

    # Parse the nested identifier (format: "inner_var:outer_var")
    id_parts <- strsplit(nested_id, ":")[[1]]
    inner_id <- id_parts[1]
    outer_id <- id_parts[2]

    # Find matching rows
    is_nested <- (fixed[[inner_var]] == inner_id &
      fixed[[outer_var]] == outer_id)

    fixed[is_nested, response_var] <-
      fixed[is_nested, response_var] - nested_intercept
  }

  fixed
}


#' Subtract random effects from data
#'
#' In order to accurately estimate the local effect size of an interaction,
#' we need to fit standard linear models that have the random effects
#' removed from the data.
#'
#' @param d The data with random effects included, used to fit model m
#' @param m The model with random effects structure that must be removed
#' @return Data frame with random effects subtracted
#' @importFrom utils tail
#' @export
subtract_random_effects <- function(
  d, # The data with random effects included, used to fit model m.
  m # The model with random effects structure that must be removed.
) {
  random_effects <- lme4::ranef(m)
  re_names <- names(random_effects)

  if (length(re_names) == 1 && !grepl(":", re_names[1])) {
    fixed <- subtract_simple_random_effects(d, m)
  } else {
    fixed <- subtract_nested_random_effects(d, m)
  }
  fixed
}

#' Strip random effects from formula
#'
#' Reconstruct a formula with only fixed effects.
#' Example A: y ~ x * z + (1|subject) -> y ~ x * z
#' Example B: y ~ x + (1|subject) -> y ~ x
#' Example C: y ~ (1|subject) -> y ~ 1
#' Example D: y ~ (1|subject) + (1|group) -> y ~ 1
#' For most formulas, this will be equivalent to nobars(formula_obj).
#'
#' @param formula_obj Formula object with random effects
#' @return Formula object with only fixed effects
#' @export
strip_random_effects <- function(formula_obj) {
  # Split formula into components
  terms <- as.character(formula_obj)

  # Get left-hand side (response variable)
  lhs <- terms[2]

  # Get right-hand side (fixed effects)
  rhs <- terms[3]

  # Remove random effects terms (anything containing " | ")
  fixed_terms <- unlist(strsplit(rhs, " \\+ "))
  fixed_terms <- fixed_terms[!grepl(" \\| ", fixed_terms)]

  # If no fixed effects remain, use "1" as the RHS (just an intercept).
  if (length(fixed_terms) == 0) {
    fixed_rhs <- "1"
  } else {
    fixed_rhs <- paste(fixed_terms, collapse = " + ")
  }

  # Reconstruct formula
  new_formula <- stats::as.formula(paste(lhs, "~", fixed_rhs))
  new_formula
}

#' Calculate Cohen's local f-squared
#'
#' For estimation statistics, Cohen's f2 (f-squared) is a measure of effect
#' size. The local version of the f2 statistic quantifies the amount of variance
#' accounted for by variables of interest, relative to the amount of variance
#' left unaccounted for by the model.
#'
#' @param model_a The model with the effect of interest included
#' @param model_b The model with the effect of interest removed
#' @return Cohen's f-squared value
#' @export
cohens_local_fsquared <- function(model_a, model_b) {
  r2_a <- summary(model_a)$r.squared
  r2_b <- summary(model_b)$r.squared
  f2 <- (r2_a - r2_b) / (1 - r2_a)
  f2
}

#' Subtract random effects and calculate f-squared
#'
#' There is a caveat to computing f2 for LME models; when we remove the effect
#' of interest, it is possible that the reduced model is estimated with a
#' drastically different random effect structure, misrepresenting the variance
#' accounted for by the fixed effects (Selya et al., 2012). To circumvent this
#' issue, we first fit the full model to estimate a random effect structure
#' and subtract the estimated random effects from the data. We then fit two
#' standard linear models to the adjusted data (one full, one reduced) and use
#' their R2 values to calculate f2. In this way, we enforce that the full and
#' reduced models have identical random effect structure, and that f2 explicitly
#' captures the variance accounted for by the fixed effects.
#'
#' @param d Data frame
#' @param model_a The model with the effect of interest included, and the source
#' of random effects to be subtracted from the data
#' @param model_b The model with the effect of interest removed
#' @return List with adjusted data, models, and f-squared value
#' @export
subtract_ranef_get_fsquared <- function(d, model_a, model_b) {
  d_fixef <- subtract_random_effects(d, model_a)

  # Remove random effects structure from the formulas, and refit using the
  # data that also has random effects structure removed.
  # Setting the model's call attribute ensures that it prints nicely.
  fa_fixef <- strip_random_effects(stats::formula(model_a))
  ma_fixef <- stats::lm(fa_fixef, data = d_fixef)
  ma_fixef$call <- call("lm", formula = fa_fixef, data = quote(d_fixef))

  fb_fixef <- strip_random_effects(stats::formula(model_b))
  mb_fixef <- stats::lm(fb_fixef, data = d_fixef)
  mb_fixef$call <- call("lm", formula = fb_fixef, data = quote(d_fixef))

  fsquared <- cohens_local_fsquared(ma_fixef, mb_fixef)

  list(
    data = d_fixef,
    models = list(a = ma_fixef, b = mb_fixef),
    fsquared = fsquared
  )
}

#' Calculate Cohen's d analogue for mixed effects models
#'
#' The ratio of mean to variance is used as a measure of effect size, with the
#' combined residual variance and random effect variance in the denominator
#' (analogous to Cohen's D but for LME models).
#'
#' @param glht General linear hypothesis test object
#' @param model Mixed effects model
#' @return Matrix of Cohen's d analogue values
#' @export
cohens_d_analogue <- function(glht, model) {
  # We are being very conservative here.
  # lme4::sigma(model) would also be fine as a measure of variance.
  # There is no standard.
  denom <- sqrt(sum(as.data.frame(lme4::VarCorr(model))$vcov))
  eff_size <- as.matrix(stats::coef(glht)) / denom
  colnames(eff_size) <- "Cohen's d analogue"
  eff_size
}
