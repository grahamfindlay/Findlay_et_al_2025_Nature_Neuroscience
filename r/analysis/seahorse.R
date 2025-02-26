# Create environment for seahorse functions
library(arrow)
library(assertthat)
library(fastDummies)
library(ggplot2)
library(here)
library(lme4)
library(multcomp)
library(rlang)

sh <- env()

sh$alpha <- 0.05

sh$load_measures <- function() {
  path <- here("data", "condition_measures.pqt")
  d <- read_parquet(path)
  exp_levels <- c("Novelty", "Locomotion", "Dual")
  d$experiment <- factor(d$experiment, levels = exp_levels)
  d$subject <- factor(d$subject)
  return(d)
}

sh$load_contrasts <- function() {
  path <- here("data", "condition_contrasts.pqt")
  d <- read_parquet(path)
  exp_levels <- c("Novelty", "Locomotion", "Dual")
  d$experiment <- factor(d$experiment, levels = exp_levels)
  d$subject <- factor(d$subject)
  return(d)
}

# If the full model is:
#   Cortical SWA ~ Condition * Experiment + (1 | Subject)
# then Cortical SWA is the response variable, Condition and Experiment are
# the primary and secondary covariates, and Subject is the random factor.
# The reduced model is therefore:
#   Cortical SWA ~ Condition + Experiment + (1 | Subject)
# The null model is:
#   Cortical SWA ~ Experiment + (1 | Subject)
sh$fit_nested_models <- function(
    d, response_var, primary_covariate, secondary_covariate, random_factor) {
  # Quote the variables in case they include spaces or special characters.
  # This also keeps them from being evaluated prematurely.
  primary_covariate_ <- paste0("`", primary_covariate, "`")
  secondary_covariate_ <- paste0("`", secondary_covariate, "`")
  response_var_ <- paste0("`", response_var, "`")
  random_factor_ <- paste0("`", random_factor, "`")
  random_intercept_ <- paste0("(1 | ", random_factor_, ")")

  # Create formula strings
  f1 <- reformulate(
    c(paste(primary_covariate_, "*", secondary_covariate_), random_intercept_),
    response = response_var_
  ) # Full model

  f2 <- reformulate(
    c(primary_covariate_, secondary_covariate_, random_intercept_),
    response = response_var_
  ) # Reduced model


  f3 <- reformulate(
    c(secondary_covariate_, random_intercept_),
    response = response_var_
  ) # Null model

  # Fit models
  m1 <- lmer(f1, data = d, REML = FALSE) # Full model
  m2 <- lmer(f2, data = d, REML = FALSE) # Reduced model
  m3 <- lmer(f3, data = d, REML = FALSE) # Null model

  return(list(full = m1, reduced = m2, null = m3))
}

# In order to accurately estimate the local effect size of an interaction,
# we need to fit standard linear models that have the random effects
# removed from the data.
sh$subtract_random_effects <- function(
  d, # The data with random effects included, used to fit model m.
  m # The model with random effects structure that must be removed.
) {
  fixed <- d # Make a copy of the data

  # Extract model components. If the model is:
  #   Cortical SWA ~ Condition * Experiment + (1 | Subject)
  # then Cortical SWA is the response variable, Condition and Experiment are
  # the primary and secondary covariates, and Subject is the random factor.
  random_effects <- ranef(m)
  random_factor <- names(random_effects)
  response_var <- all.vars(formula(m))[1]
  assert_that(random_factor == tail(all.vars(formula(m)), n = 1))

  # For each random factor level (e.g., for each subject),
  # subtract the random effect from the response variable.
  random_factor_levels <- levels(fixed[[random_factor]])
  for (fctr in random_factor_levels) {
    intercept <- random_effects[[random_factor]][fctr, 1]
    is_fctr <- (fixed[[random_factor]] == fctr)
    fixed[is_fctr, response_var] <- fixed[is_fctr, response_var] - intercept
  }
  return(fixed)
}

# Reconstruct a formula with only fixed effects.
# Example A: y ~ x * z + (1|subject) -> y ~ x * z
# Example B: y ~ x + (1|subject) -> y ~ x
# Example C: y ~ (1|subject) -> y ~ 1
# Example D: y ~ (1|subject) + (1|group) -> y ~ 1
# For most formulas, this will be equivalent to nobars(formula_obj).
sh$strip_random_effects <- function(formula_obj) {
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
  new_formula <- as.formula(paste(lhs, "~", fixed_rhs))
  return(new_formula)
}

# For estimation statistics, Cohen's f2 (f-squared) is a measure of effect size.
# The local version of the f2 statistic quantifies the amount of variance
# accounted for by variables of interest, relative to the amount of variance
# left unaccounted for by the model.
#
# model_a: The model with the effect of interest included.
# model_b: The model with the effect of interest removed.
#
# For example, when estimating the effect size of an interaction,
# model_a and model_b might be the "fixed" full and reduced models
# from sh$fit_models_for_interaction_effect_size.
sh$cohens_local_fsquared <- function(model_a, model_b) {
  r2_a <- summary(model_a)$r.squared
  r2_b <- summary(model_b)$r.squared
  f2 <- (r2_a - r2_b) / (1 - r2_a)
  return(f2)
}

# There is a caveat to computing f2 for LME models; when we remove the effect
# of interest, it is possible that the reduced model is estimated with a
# drastically different random effect structure, misrepresenting the variance
# accounted for by the fixed effects (Selya et al., 2012). To circumvent this
# issue, we first fit the full model to estimate a random effect structure
# and subtract the estimated random effects from the data. We then fit two
# standard linear models to the adjusted data (one full, one reduced) and use
# their R2 values to calculate f2. In this way, we enforce that the full and
# reduced models have identical random effect structure, and that f2 explicitly
# captures the variance accounted for by the fixed effects.
#
# model_a: The model with the effect of interest included, and the source of
#   random effects to be subtracted from the data.
# model_b: The model with the effect of interest removed.
sh$subtract_ranef_and_get_fsquared <- function(d, model_a, model_b) {
  d_fixef <- sh$subtract_random_effects(d, model_a)

  # Remove random effects structure from the formulas, and refit using the
  # data that also has random effects structure removed.
  # Setting the model's call attribute ensures that it prints nicely.
  fa_fixef <- sh$strip_random_effects(formula(model_a))
  ma_fixef <- lm(fa_fixef, data = d_fixef)
  ma_fixef$call <- call("lm", formula = fa_fixef, data = quote(d_fixef))

  fb_fixef <- sh$strip_random_effects(formula(model_b))
  mb_fixef <- lm(fb_fixef, data = d_fixef)
  mb_fixef$call <- call("lm", formula = fb_fixef, data = quote(d_fixef))

  fsquared <- sh$cohens_local_fsquared(ma_fixef, mb_fixef)

  return(
    list(
      data = d_fixef,
      models = list(a = ma_fixef, b = mb_fixef),
      fsquared = fsquared
    )
  )
}

# The ratio of mean to variance is used as a measure of effect size, with the
# combined residual variance and random effect variance in the denominator
# (analogous to Cohen’s D but for LME models).
sh$cohens_d_analogue <- function(glht, model) {
  # We are being very conservative here.
  # sigma(model) would also be fine as a measure of variance.
  # There is no standard.
  denom <- sqrt(sum(as.data.frame(VarCorr(model))$vcov))
  eff_size <- as.matrix(coef(glht)) / denom
  colnames(eff_size) <- "Cohen's d analogue"
  return(eff_size)
}

sh$format_posthoc_summary <- function(posthoc) {
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
  return(msg)
}

sh$get_anova_pval <- function(anova_result) {
  return(anova_result["Pr(>Chisq)"][[1]][2])
}

# Tests for an interaction between primary and secondary covariates.
# In the event of a significant interaction, this function performs
# optional post-hoc tests specified by the single contrast matrix provided.
# If you need multiple contrast matrices (e.g., one per experiment),
# use sh$test_interaction_with_experiment() below.
sh$test_interaction <- function(dat, models, contrast_matrix = NULL) {
  # Test for interaction
  interaction <- list()
  interaction$anova <- anova(models$full, models$reduced)
  interaction$pval <- sh$get_anova_pval(interaction$anova)

  # If interaction is significant, estimate its effect size
  # and perform post-hoc tests
  if (interaction$pval < sh$alpha) {
    # Estimate interaction effect size
    interaction$effect_size <-
      sh$subtract_ranef_and_get_fsquared(dat, models$full, models$reduced)

    # Perform optional post-hoc tests for nonzero slopes
    # Also get confidence intervals and effect sizes
    if (!is.null(contrast_matrix)) {
      ph <- list()
      ph$glht <- glht(models$full, contrast_matrix)
      ph$ci <- confint(ph$glht)
      ph$effect_size <- sh$cohens_d_analogue(ph$glht, models$full)
      interaction$posthoc <- ph
    }
  }
  return(interaction)
}

# If the main effect is significant and a contrast matrix is provided,
# this function performs optional post hoc tests.
# This function does nothing to stop you from testing for a main effect
# when such a test is incoherent (e.g., when a significant interaction
# exists). That's on you.
sh$test_main_effect <- function(dat, models, contrast_matrix = NULL) {
  # Test for main effect
  main_effect <- list()
  main_effect$anova <- anova(models$reduced, models$null)
  main_effect$pval <- sh$get_anova_pval(main_effect$anova)

  # If main effect is significant, estimate its effect size
  # and perform post-hoc tests
  if (main_effect$pval < sh$alpha) {
    # Estimate main effect size
    main_effect$effect_size <- sh$subtract_ranef_and_get_fsquared(
      dat,
      models$reduced,
      models$null
    )

    # Perform optional post-hoc tests
    if (!is.null(contrast_matrix)) {
      ph <- list()
      ph$glht <- glht(models$reduced, contrast_matrix)
      ph$ci <- confint(ph$glht)
      ph$effect_size <- sh$cohens_d_analogue(ph$glht, models$reduced)
      # Attach post-hoc tests to main effect
      main_effect$posthoc <- ph
    }
  }
  return(main_effect)
}

# Tests for an interaction between any categorical primary covariate
# (e.g., condition, or state) and experiment (as a secondary covariate).
# You can also test for interactions with experiment using the more generic
# sh$test_interaction function above, but this function is more convenient
# for the common case where you have a contrast matrix for each experiment.
sh$test_interaction_with_experiment <- function(
  dat,
  models,
  novelty_matrix,
  locomotion_matrix,
  dual_matrix
) {
  # Test for interaction
  interaction <- list()
  interaction$anova <- anova(models$full, models$reduced)
  interaction$pval <- sh$get_anova_pval(interaction$anova)

  # If interaction is significant, estimate its effect size
  # and perform post-hoc tests
  if (interaction$pval < sh$alpha) {
    # Estimate interaction effect size
    interaction$effect_size <- sh$subtract_ranef_and_get_fsquared(
      dat,
      models$full,
      models$reduced
    )

    # Perform post-hoc tests
    ph <- list()
    # Novelty post-hoc tests, confidence intervals, and effect sizes
    ph$nod$glht <- glht(models$full, novelty_matrix)
    ph$nod$ci <- confint(ph$nod$glht)
    ph$nod$effect_size <- sh$cohens_d_analogue(ph$nod$glht, models$full)
    # Locomotion post-hoc tests, confidence intervals, and effect sizes
    ph$cow$glht <- glht(models$full, locomotion_matrix)
    ph$cow$ci <- confint(ph$cow$glht)
    ph$cow$effect_size <- sh$cohens_d_analogue(ph$cow$glht, models$full)
    # Dual post-hoc tests, confidence intervals, and effect sizes
    ph$ctn$glht <- glht(models$full, dual_matrix)
    ph$ctn$ci <- confint(ph$ctn$glht)
    ph$ctn$effect_size <- sh$cohens_d_analogue(ph$ctn$glht, models$full)
    # Attach post-hoc tests to interaction
    interaction$posthoc <- ph
  }
  return(interaction)
}

# Homeostasis testing ----------------------------------------------------------

sh$load_measures_for_homeostasis_testing <- function() {
  d <- sh$load_measures()

  condition_levels <- c(
    "Early BSL NREM",
    "Late BSL NREM",
    "Early EXT Wake",
    "Late EXT Wake",
    "Early REC NREM",
    "Late REC NREM"
  )
  d <- subset(d, d$condition %in% condition_levels)
  d$condition <- factor(d$condition, levels = condition_levels)

  return(d)
}

sh$homeostasis_novelty_matrix <- rbind(
  "E_REC - L_BSL" = c(0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "L_REC - E_REC" = c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "L_EWK - E_EWK" = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "E_BSL - L_BSL" = c(0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "E_REC - E_BSL" = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)

sh$homeostasis_locomotion_matrix <- rbind(
  "E_REC - L_BSL" = c(0, -1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0),
  "L_REC - E_REC" = c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0),
  "L_EWK - E_EWK" = c(0, 0, -1, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0),
  "E_BSL - L_BSL" = c(0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "E_REC - E_BSL" = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
)

sh$homeostasis_dual_matrix <- rbind(
  "E_REC - L_BSL" = c(0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0),
  "L_REC - E_REC" = c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1),
  "L_EWK - E_EWK" = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0),
  "E_BSL - L_BSL" = c(0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0),
  "E_REC - E_BSL" = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
)

sh$homeostasis_main_effect_matrix <- rbind(
  "E_REC - L_BSL" = c(0, -1, 0, 0, 1, 0, 0, 0),
  "L_REC - E_REC" = c(0, 0, 0, 0, -1, 1, 0, 0),
  "L_EWK - E_EWK" = c(0, 0, -1, 1, 0, 0, 0, 0),
  "E_BSL - L_BSL" = c(0, -1, 0, 0, 0, 0, 0, 0),
  "E_REC - E_BSL" = c(0, 0, 0, 0, 1, 0, 0, 0)
)

sh$test_homeostasis <- function(response_variable) {
  r <- list()
  r$data <- sh$load_measures_for_homeostasis_testing()
  r$models <- sh$fit_nested_models(
    r$data,
    response_variable,
    "condition",
    "experiment",
    "subject"
  )

  r$interaction <- sh$test_interaction_with_experiment(
    r$data,
    r$models,
    sh$homeostasis_novelty_matrix,
    sh$homeostasis_locomotion_matrix,
    sh$homeostasis_dual_matrix
  )
  if (r$interaction$pval >= sh$alpha) {
    r$main_effect <- sh$test_main_effect(
      r$data,
      r$models,
      sh$homeostasis_main_effect_matrix
    )
  }
  return(r)
}

# Correlation testing ----------------------------------------------------------

# If x_condition = "ExtWake", x_variable = "Theta",
# y_condition = "RecNREM", y_variable = "SPWs"...
#
# Transform data d of the form:
#
# subject experiment condition Theta SPWs
# 1       A          ExtWake   10    20
# 1       A          RecNREM   15    25
# 1       B          ExtWake   20    30
# 1       B          RecNREM   25    35
#
# into dxy of the form:
#
# subject experiment `ExtWake Theta` `RecNREM SPWs`
# 1       A          10              25
# 1       B          20              35
# 2       A          ...             ...
#
# where xvar = "`ExtWake Theta`" and yvar = "`RecNREM SPWs`"
sh$get_correlation_data <- function(
  d,
  x_condition,
  x_variable,
  y_condition,
  y_variable
) {
  xvar <- paste(x_condition, x_variable)
  dx <- d[d$condition == x_condition, ]
  dx[xvar] <- dx[x_variable]
  dx <- dx[, c("subject", "experiment", xvar)]

  yvar <- paste(y_condition, y_variable)
  dy <- d[d$condition == y_condition, ]
  dy[yvar] <- dy[y_variable]
  dy <- dy[, c("subject", "experiment", yvar)]

  dxy <- merge(dx, dy, by = c("subject", "experiment"))
  return(list(data = dxy, xvar = xvar, yvar = yvar))
}

# If you have covariate X (hippocampal theta-delta ratio),
# categorical variable Z (experiment) with levels A (NOD), B (COW), C (CTN).
#
# The full model y ~ X * Z will have terms [intercept, zB, zC, X, zB:X, zC:X]
# To test for an effect for A (Novelty), you test if X = 0
# To test for an effect for B (Locomotion), you test if X + zB:X = 0

# To test for an effect for C (Locomotion + Novelty), you test if X + zC:X = 0
sh$correlation_interaction_matrix <- rbind(
  "Novelty" = c(0, 1, 0, 0, 0, 0),
  "Locomotion" = c(0, 1, 0, 0, 1, 0),
  "Dual" = c(0, 1, 0, 0, 0, 1)
)

sh$fit_dummy_models_for_fsquared <- function(
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
  ) # TODO: These should probably be quoted for safety.

  # Quote the variables in case they include spaces or special characters.
  # This also keeps them from being evaluated prematurely.
  primary_covariate_ <- paste0("`", primary_covariate, "`")
  response_var_ <- paste0("`", response_var, "`")
  random_factor_ <- paste0("`", random_factor, "`")
  random_intercept_ <- paste0("(1 | ", random_factor_, ")")

  # Create formula string for full dummy model
  # The "0" forces the intercept to be zero, which allows the dummy variables
  # to appear in the model, which would otherwise be overdetermined.
  simple_fixefs <- paste(dummy_vars, collapse = " + ")
  intractn_terms <- paste0(dummy_vars, ":", primary_covariate_)
  intractn_fixefs <- paste(intractn_terms, collapse = " + ")
  dummy_formula <- reformulate(
    c(
      paste("0", "+", simple_fixefs, "+", intractn_fixefs),
      random_intercept_
    ),
    response = response_var_
  )

  # Fit full dummy model
  dummy_model <- lmer(dummy_formula, data = dummy_d, REML = FALSE)

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
    reduced_formula <- reformulate(
      c(
        paste("0", "+", simple_fixefs, "+", remaining_fixefs),
        random_intercept_
      ),
      response = response_var_
    )

    # Fit reduced model and store it
    reduced_models[[dummy]] <- lmer(
      reduced_formula,
      data = dummy_d,
      REML = FALSE
    )
  }

  return(
    list(
      data = dummy_d,
      models = list(
        full = dummy_model,
        reduced = reduced_models
      )
    )
  )
}

# Get Cohen's f-squared for posthoc tests, when Cohen's d is not appropriate.
# The issue is that Cohen's d estimates may depend on the scale of the inputs.
# You can change the scale of the inputs, making the estimates larger or
# smaller arbitrarily, and the residual error won't change (linear
# transformation of a linear model). So that notion of effect size doesn't
# work, but something like Cohen's f-squared does.
sh$get_fsquared_for_interaction_posthocs <- function(
  d,
  response_var,
  primary_covariate,
  secondary_covariate,
  random_factor
) {
  dummy <- sh$fit_dummy_models_for_fsquared(
    d,
    response_var,
    primary_covariate,
    secondary_covariate,
    random_factor
  )

  effect_sizes <- list()
  for (lvl in names(dummy$models$reduced)) {
    effect_sizes[[lvl]] <- sh$subtract_ranef_and_get_fsquared(
      dummy$data,
      dummy$models$full,
      dummy$models$reduced[[lvl]]
    )
  }
  return(effect_sizes)
}

sh$test_condition_correlation <- function(
  x_condition,
  x_variable,
  y_condition,
  y_variable
) {
  measures <- sh$load_measures()
  dxy <- sh$get_correlation_data(
    measures,
    x_condition,
    x_variable,
    y_condition,
    y_variable
  )

  r <- list() # Hold results
  r$data <- dxy$data
  r$models <- sh$fit_nested_models(
    dxy$data,
    dxy$yvar,
    dxy$xvar,
    "experiment", 
    "subject"
  )

  r$interaction <- sh$test_interaction(
    r$data,
    r$models,
    sh$correlation_interaction_matrix
  )
  r$interaction$posthoc$effect_size <- sh$get_fsquared_for_interaction_posthocs(
    r$data,
    dxy$yvar,
    dxy$xvar,
    "experiment",
    "subject"
  ) # Overwrite faulty Cohen's d analogues with f-squared
  assert_that(
    all(
      rownames(sh$correlation_interaction_matrix)
      == names(r$interaction$posthoc$effect_size)
    ),
    msg = "Contrast matrix rows and posthoc levels do not exactly match!"
  ) # Not strictly necessary, but nice for later formatting of outputs.
  if (r$interaction$pval >= sh$alpha) {
    r$main_effect <- sh$test_main_effect(r$data, r$models)
  }
  return(r)
}

sh$test_contrast_correlation <- function(x_variable, y_variable) {
  r <- list()
  r$data <- sh$load_contrasts()
  r$models <- sh$fit_nested_models(
    r$data,
    y_variable,
    x_variable,
    "experiment",
    "subject"
  )

  r$interaction <- sh$test_interaction(
    r$data,
    r$models,
    sh$correlation_interaction_matrix
  )
  r$interaction$posthoc$effect_size <- sh$get_fsquared_for_interaction_posthocs(
    r$data,
    y_variable,
    x_variable,
    "experiment",
    "subject"
  ) # Overwrite faulty Cohen's d analogues with f-squared
  assert_that(
    all(
      rownames(sh$correlation_interaction_matrix)
      == names(r$interaction$posthoc$effect_size)
    ),
    msg = "Contrast matrix rows and posthoc levels do not exactly match!"
  ) # Not strictly necessary, but nice for later formatting of outputs.
  if (r$interaction$pval >= sh$alpha) {
    r$main_effect <- sh$test_main_effect(r$data, r$models)
  }
  return(r)
}

# SG/FG testing ---------------------------------------------------------------

sh$load_sgfg_by_epoch_type <- function() {
  path <- here("data", "AeryJones_sgfg_medians_by_epoch_type.pqt")
  d <- read_parquet(path)

  # Ensure proper level ordering, for contrast matrix coding
  exp_levels <- c("Novelty", "Locomotion", "Dual")
  d$Experiment <- factor(d$Experiment, levels = exp_levels)

  epoch_levels <- c("SPW Wake", "No-SPW Wake", "SPW NREM", "No-SPW NREM")
  d$`Epoch Type` <- factor(d$`Epoch Type`, levels = epoch_levels)

  d$Subject <- factor(d$Subject)

  d <- d[d$ROI == "CA1-slm", ]
  return(d)
}

# Extended Wake change testing -------------------------------------------------

sh$load_measures_for_ewk_change_testing <- function() {
  d <- sh$load_measures()

  condition_levels <- c("Early EXT Wake", "Late EXT Wake")
  d <- subset(d, d$condition %in% condition_levels)
  d$condition <- factor(d$condition, levels = condition_levels)

  return(d)
}

sh$ewk_change_novelty_matrix <-
  rbind("L_EWK - E_EWK" = c(0, 1, 0, 0, 0, 0))
sh$ewk_change_locomotion_matrix <-
  rbind("L_EWK - E_EWK" = c(0, 1, 0, 1, 0, 0))
sh$ewk_change_dual_matrix <-
  rbind("L_EWK - E_EWK" = c(0, 1, 0, 0, 0, 1))
sh$ewk_change_main_effect_matrix <-
  rbind("L_EWK - E_EWK" = c(0, 1))

sh$test_ewk_change <- function(response_variable) {
  r <- list()
  r$data <- sh$load_measures_for_ewk_change_testing()
  r$models <- sh$fit_nested_models(
    r$data,
    response_variable,
    "condition",
    "experiment",
    "subject"
  )

  r$interaction <- sh$test_interaction_with_experiment(
    r$data,
    r$models,
    sh$ewk_change_novelty_matrix,
    sh$ewk_change_locomotion_matrix,
    sh$ewk_change_dual_matrix
  )
  if (r$interaction$pval >= sh$alpha) {
    r$main_effect <- sh$test_main_effect(
      r$data,
      r$models,
      sh$ewk_change_main_effect_matrix
    )
  }
  return(r)
}

# Novelty vs Dual first exposure testing ---------------------------------------
sh$load_measures_for_first_novelty_testing <- function(only_dual_subjects) {
  d <- sh$load_measures()
  d <- droplevels(subset(d, d$condition %in% c("Early NOD Wake")))
  d <- droplevels(subset(d, d$experiment %in% c("Novelty", "Dual")))
  if (only_dual_subjects) {
    dual_subjects <- c(
      "CNPIX15-Claude",
      "CNPIX17-Hans",
      "CNPIX18-Pier",
      "CNPIX19-Otto",
      "CNPIX20-Ernst"
    )
    d <- droplevels(subset(d, d$subject %in% dual_subjects))
  }
  return(d)
}

sh$fit_basic_models <- function(
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
  f1 <- reformulate(
    c(primary_covariate, random_intercept),
    response = response_var
  ) # Full model

  f2 <- reformulate(
    c(random_intercept),
    response = response_var
  ) # Null model


  # Fit models
  m1 <- lmer(f1, data = d, REML = FALSE) # Full model
  m2 <- lmer(f2, data = d, REML = FALSE) # Null model

  return(list(full = m1, null = m2))
}

# Ripple frequency testing -----------------------------------------------------

sh$load_ripple_frequency_data <- function() {
  path <- here("data", "ripple_frequency_by_state.pqt")
  d <- read_parquet(path)

  exp_levels <- c("Novelty", "Locomotion", "Dual")
  d$experiment <- factor(d$experiment, levels = exp_levels)

  state_levels <- c("Wake", "NREM")
  d$state <- factor(d$state, state_levels)

  d$subject <- factor(d$subject)

  return(d)
}

# Occupancy testing ------------------------------------------------------------

sh$load_occupancy_data <- function() {
  # Eugene has already been excluded.
  path <- here("data", "sleep_period_fractional_occupancy.pqt")
  d <- read_parquet(path)

  exp_levels <- c("Novelty", "Locomotion", "Dual")
  d$experiment <- factor(d$experiment, levels = exp_levels)

  day_levels <- c("Baseline", "Recovery")
  d$day <- factor(d$day, levels = day_levels)

  state_levels <- c("NREM", "IS", "REM", "Wake", "MA")
  d$state <- factor(d$state, levels = state_levels)

  d$subject <- factor(d$subject)

  d <- droplevels(subset(d, d$experiment %in% c("Novelty", "Locomotion")))
  return(d)
}

sh$single_day_occupancy_contrast_matrix <- rbind(
  "NREM: COW - NOD" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  "IS: COW - NOD"   = c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
  "REM: COW - NOD"  = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0),
  "Wake: COW - NOD" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0),
  "MA: COW - NOD"   = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1)
)

sh$single_exp_occupancy_contrast_matrix <- rbind(
  "NREM: REC - BSL" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  "IS: REC - BSL"   = c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
  "REM: REC - BSL"  = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0),
  "Wake: REC - BSL" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0),
  "MA: REC - BSL"   = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1)
)

sh$test_single_day_occupancy <- function(day) {
  d <- sh$load_occupancy_data()

  r <- list()
  r$data <- droplevels(d[d$day == day, ])
  r$models <- sh$fit_nested_models(
    r$data,
    "fractional_occupancy",
    "experiment",
    "state",
    "subject"
  )

  r$interaction <- sh$test_interaction(
    r$data,
    r$models,
    sh$single_day_occupancy_contrast_matrix
  )
  if (r$interaction$pval >= sh$alpha) {
    r$main_effect <- sh$test_main_effect(r$data, r$models)
    assert_that(
      r$main_effect$pval >= sh$alpha,
      msg = "Contrast matrix for main effect posthocs not implemented!"
    )
  }
  return(r)
}

sh$test_single_experiment_occupancy <- function(experiment) {
  d <- sh$load_occupancy_data()

  r <- list()
  r$data <- droplevels(d[d$experiment == experiment, ])
  r$models <- sh$fit_nested_models(
    r$data,
    "fractional_occupancy",
    "day",
    "state",
    "subject"
  )

  r$interaction <- sh$test_interaction(
    r$data,
    r$models,
    sh$single_exp_occupancy_contrast_matrix
  )
  if (r$interaction$pval >= sh$alpha) {
    r$main_effect <- sh$test_main_effect(r$data, r$models)
    assert_that(
      r$main_effect$pval >= sh$alpha,
      msg = "Contrast matrix for main effect posthocs not implemented!"
    )
  }
  return(r)
}

sh$load_rec_rem_vs_ewk_tdr_data <- function() {
  path <- here("data", "rec_rem_vs_ewk_tdr.pqt")
  dat <- read_parquet(path)

  dat <- dat[complete.cases(dat), ] # Drop (intentionally) missing values

  exp_levels <- c("Novelty", "Locomotion", "Dual")
  dat$Experiment <- factor(dat$Experiment, levels = exp_levels)

  dat$Subject <- factor(dat$Subject)

  return(dat)
}
