#' Create contrast matrix for two-way interaction simple slopes
#'
#' This function generates a contrast matrix for testing simple slopes of a
#' primary variable at different levels of a secondary variable in the context
#' of a two-way interaction.
#'
#' @param model A fitted model object (typically from lme4) containing the
#'   interaction terms
#' @param primary_var Character string specifying the name of the primary
#'   variable (the variable whose simple slopes are being tested)
#' @param secondary_var Character string specifying the name of the secondary
#'   variable (the moderator variable)
#'
#' @return A contrast matrix (k) with rows corresponding to levels of the
#'   secondary variable and columns corresponding to model coefficients. Each
#'   row represents the linear combination of coefficients needed to test the
#'   simple slope of the primary variable at that level of the secondary
#'   variable.
make_two_way_contrasts <- function(
  model,
  primary_var,
  secondary_var
) {
  coef_names <- names(lme4::fixef(model))
  base_slope_name <- primary_var

  # Only as they appear in the model, in case any that should be there are
  # missing due to rank-deficiency.
  interaction_terms <- grep(
    paste0("^", primary_var, ":", secondary_var),
    coef_names,
    value = TRUE
  )
  interaction_levels <- sub(
    paste0("^", primary_var, ":", secondary_var),
    "",
    interaction_terms
  )

  # Only as they appear in the model, in case any that should be there are
  # missing due to rank-deficiency.
  cov_terms <- grep(paste0("^", secondary_var), coef_names, value = TRUE)
  cov_levels <- sub(paste0("^", secondary_var), "", cov_terms)

  # As they appear in the data.
  all_cov_levels <- levels(model@frame[[secondary_var]])
  ref_level <- setdiff(all_cov_levels, cov_levels)

  # Append reference level to interaction terms to ensure it's included
  interaction_levels <- c(ref_level, interaction_levels)

  k <- matrix(
    0,
    nrow = length(interaction_levels),
    ncol = length(coef_names)
  )
  rownames(k) <- interaction_levels
  colnames(k) <- coef_names
  for (i in seq_along(interaction_levels)) {
    level <- interaction_levels[i]
    k[i, base_slope_name] <- 1
    if (level != ref_level) {
      int_name <- paste0(primary_var, ":", secondary_var, level)
      if (int_name %in% interaction_terms) {
        k[i, int_name] <- 1
      } else {
        stop(
          paste("Interaction term", int_name, "not found in model coefficients")
        )
      }
    }
  }
  k
}

#' Create contrast matrix for three-way interaction simple slopes
#'
#' This function generates a contrast matrix for testing simple slopes of a
#' primary variable at different combinations of levels of secondary and
#' tertiary variables in the context of a three-way interaction.
#'
#' @param model A fitted model object (typically from lme4) containing the
#'   three-way interaction terms
#' @param primary_var Character string specifying the name of the primary
#'   variable (the variable whose simple slopes are being tested)
#' @param secondary_var Character string specifying the name of the secondary
#'   moderator variable
#' @param tertiary_var Character string specifying the name of the tertiary
#'   moderator variable
#'
#' @return A contrast matrix (k) with rows corresponding to combinations of
#'   levels of the secondary and tertiary variables, and columns corresponding
#'   to model coefficients. Each row represents the linear combination of
#'   coefficients needed to test the simple slope of the primary variable at
#'   that combination of moderator levels. Rows with missing interaction terms
#'   are excluded with warnings.
make_three_way_contrasts <- function(
  model,
  primary_var,
  secondary_var,
  tertiary_var
) {
  # Extract coefficient names once
  coef_names <- names(lme4::fixef(model))
  base_slope_name <- primary_var

  # Helpers --------------------------------------------------------------------
  # Non-reference levels for a covariate, restricting to *main-effect* terms
  get_nonref_levels <- function(var) {
    terms <- grep(paste0("^", var), coef_names, value = TRUE)
    # Keep only main-effect terms (no colons)
    terms <- terms[!grepl(":", terms, fixed = TRUE)]
    sub(paste0("^", var), "", terms)
  }

  # Present levels in the data (drop empty ones)
  sec_levels_all <- levels(droplevels(model@frame[[secondary_var]]))
  ter_levels_all <- levels(droplevels(model@frame[[tertiary_var]]))

  # Determine reference levels -------------------------------------------------
  sec_nonref <- get_nonref_levels(secondary_var)
  ter_nonref <- get_nonref_levels(tertiary_var)

  sec_ref <- setdiff(sec_levels_all, sec_nonref)
  ter_ref <- setdiff(ter_levels_all, ter_nonref)

  assertthat::assert_that(
    length(sec_ref) == 1 && length(ter_ref) == 1,
    msg = "Unable to unambiguously determine reference levels"
  )

  sec_ref <- sec_ref[1]
  ter_ref <- ter_ref[1]

  # Identify interaction coefficient names ------------------------------------
  # 2-way interactions that involve the primary + secondary (exclude tertiary)
  sec_int_terms <- grep(
    paste0("^", primary_var, ":", secondary_var),
    coef_names,
    value = TRUE
  )
  sec_int_terms <- sec_int_terms[
    !grepl(paste0(":", tertiary_var), sec_int_terms)
  ]

  # 2-way interactions that involve the primary + tertiary (exclude secondary)
  ter_int_terms <- grep(
    paste0("^", primary_var, ":", tertiary_var),
    coef_names,
    value = TRUE
  )
  ter_int_terms <- ter_int_terms[
    !grepl(paste0(":", secondary_var), ter_int_terms)
  ]

  # All 3-way interactions that contain primary + secondary + tertiary
  triple_int_terms <- grep(
    paste0(
      "^",
      primary_var,
      ":.*",
      secondary_var,
      ".*:.*",
      tertiary_var,
      "|^",
      primary_var,
      ":.*",
      tertiary_var,
      ".*:.*",
      secondary_var
    ),
    coef_names,
    value = TRUE,
    perl = TRUE
  )

  # Contrast matrix skeleton ---------------------------------------------------
  row_names <- as.vector(outer(
    sec_levels_all,
    ter_levels_all,
    FUN = function(a, b) paste0(a, ":", b)
  ))
  k <- matrix(0, nrow = length(row_names), ncol = length(coef_names))
  rownames(k) <- row_names
  colnames(k) <- coef_names

  keep <- logical(length(row_names))
  keep[] <- TRUE
  names(keep) <- row_names
  for (r in row_names) {
    # Base slope always gets a 1
    k[r, base_slope_name] <- 1

    # Parse the row name to extract secondary and tertiary levels
    parts <- strsplit(r, ":", fixed = TRUE)[[1]]
    s <- parts[1]
    t <- parts[2]

    # Add secondary 2-way interaction if needed
    if (s != sec_ref) {
      term <- paste0(primary_var, ":", secondary_var, s)
      if (term %in% sec_int_terms) {
        k[r, term] <- 1
      } else {
        warning(
          paste(
            "2-way interaction term",
            term,
            "not found in model coefficients"
          )
        )
        keep[r] <- FALSE
      }
    }

    # Add tertiary 2-way interaction if needed
    if (t != ter_ref) {
      term <- paste0(primary_var, ":", tertiary_var, t)
      if (term %in% ter_int_terms) {
        k[r, term] <- 1
      } else {
        warning(
          paste(
            "2-way interaction term",
            term,
            "not found in model coefficients"
          )
        )
        keep[r] <- FALSE
      }
    }

    # Add 3-way interaction if both modifiers are non-reference
    if (s != sec_ref && t != ter_ref) {
      cand1 <- paste0(primary_var, ":", secondary_var, s, ":", tertiary_var, t)
      cand2 <- paste0(primary_var, ":", tertiary_var, t, ":", secondary_var, s)
      term <- intersect(c(cand1, cand2), triple_int_terms)
      if (length(term) == 1) {
        k[r, term] <- 1
      } else {
        warning(
          paste(
            "3-way interaction term",
            cand1,
            "or",
            cand2,
            "not found in model coefficients"
          )
        )
        keep[r] <- FALSE
      }
    }
  }
  k <- k[keep, , drop = FALSE]
}

#' Test individual two-way interactions in a three-way interaction model
#'
#' This function tests the significance of individual two-way interactions by
#' comparing a full three-way interaction model against reduced models that
#' exclude specific two-way interaction terms. It performs likelihood ratio
#' tests for three comparisons: removing primary-tertiary, secondary-tertiary,
#' and primary-secondary interactions.
#'
#' @param full_minus_pst A fitted model object containing the full three-way
#'   interaction (but already excluding the three-way interaction term)
#' @param primary_covariate_ Character string specifying the primary covariate
#'   name
#' @param secondary_covariate_ Character string specifying the secondary
#'   covariate name
#' @param tertiary_covariate_ Character string specifying the tertiary
#'   covariate name
#' @param response_var_ Character string specifying the response variable name
#' @param random_intercept_ Character string specifying the random intercept
#'   term (e.g., "(1|subject)")
#'
#' @return A list containing three elements (minus_pt, minus_st, minus_ps),
#'   each containing the reduced model formula, fitted model object, ANOVA
#'   comparison, and p-value for testing the significance of the removed
#'   two-way interaction.
#'   minus_pt = primary:tertiary removed.
#'   minus_st = secondary:tertiary removed.
#'   minus_ps = primary:secondary removed.
test_two_ways <- function(
  full_minus_pst,
  primary_covariate_,
  secondary_covariate_,
  tertiary_covariate_,
  response_var_,
  random_intercept_
) {
  minus_pt_formula <- stats::reformulate(
    c(
      paste(
        primary_covariate_,
        "*",
        secondary_covariate_,
        "*",
        tertiary_covariate_,
        "-",
        primary_covariate_,
        ":",
        secondary_covariate_,
        ":",
        tertiary_covariate_,
        "-",
        primary_covariate_,
        ":",
        tertiary_covariate_
      ),
      random_intercept_
    ),
    response = response_var_
  )

  minus_st_formula <- stats::reformulate(
    c(
      paste(
        primary_covariate_,
        "*",
        secondary_covariate_,
        "*",
        tertiary_covariate_,
        "-",
        primary_covariate_,
        ":",
        secondary_covariate_,
        ":",
        tertiary_covariate_,
        "-",
        secondary_covariate_,
        ":",
        tertiary_covariate_
      ),
      random_intercept_
    ),
    response = response_var_
  )

  minus_ps_formula <- stats::reformulate(
    c(
      paste(
        primary_covariate_,
        "*",
        secondary_covariate_,
        "*",
        tertiary_covariate_,
        "-",
        primary_covariate_,
        ":",
        secondary_covariate_,
        ":",
        tertiary_covariate_,
        "-",
        primary_covariate_,
        ":",
        secondary_covariate_
      ),
      random_intercept_
    ),
    response = response_var_
  )

  # This is an utterly insane hack to avoid bizzare behavior in lme4/anova.
  # anova() will complain if the two models were not fit to the same data.
  # Makes sense. Except...
  # It checks this using the name of the variable used during model fitting.
  # Not the values.
  # https://github.com/lme4/lme4/issues/622
  r <- list()
  r$data <- full_minus_pst@frame

  minus_pt_model <- lme4::lmer(minus_pt_formula, data = r$data, REML = FALSE)
  minus_st_model <- lme4::lmer(minus_st_formula, data = r$data, REML = FALSE)
  minus_ps_model <- lme4::lmer(minus_ps_formula, data = r$data, REML = FALSE)

  list(
    minus_pt = test_interaction(
      r$data,
      list(full = full_minus_pst, reduced = minus_pt_model)
    ),
    minus_st = test_interaction(
      r$data,
      list(full = full_minus_pst, reduced = minus_st_model)
    ),
    minus_ps = test_interaction(
      r$data,
      list(full = full_minus_pst, reduced = minus_ps_model)
    )
  )
}
