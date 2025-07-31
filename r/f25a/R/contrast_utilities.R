#' Contrast utility functions
#'
#' Functions for creating custom contrast matrices for specific comparisons.

#' Create contrast matrices for pairwise comparisons of a primary covariate
#' at a given level of a secondary covariate.
#' For example, if the model is y ~ state * experiment + (1 | subject),
#' this function will create contrast matrices for testing all pairwise
#' comparisons between states within a specific experiment level
#' (e.g., "Novelty").
#'
#' @param model A fitted model object containing primary_var*secondary_var
#'   interaction.
#' @param primary_var Character string specifying the primary variable name
#' @param secondary_var Character string specifying the secondary variable name
#' @param secondary_level Character string specifying which secondary level
#'   to create contrasts for (e.g., "Novelty")
#' @return A contrast matrix with rows for each pairwise primary covariate
#'   and columns for model coefficients
#' @export
create_pairwise_primary_contrasts <- function(
  model,
  primary_var,
  secondary_var,
  secondary_level
) {
  # Get model coefficients and check inputs
  coef_names <- names(lme4::fixef(model))

  # Get primary and secondary levels from the model data
  primary_levels <- levels(model@frame[[primary_var]])
  secondary_levels <- levels(model@frame[[secondary_var]])

  # Validate secondary level
  if (!secondary_level %in% secondary_levels) {
    stop(paste(
      "secondary_level must be one of:",
      paste(secondary_levels, collapse = ", ")
    ))
  }

  # Create all pairwise combinations of primary levels
  primary_pairs <- utils::combn(primary_levels, 2, simplify = FALSE)
  n_contrasts <- length(primary_pairs)
  n_coefs <- length(coef_names)

  # Initialize contrast matrix
  contrast_matrix <- matrix(
    0,
    nrow = n_contrasts,
    ncol = n_coefs
  )
  colnames(contrast_matrix) <- coef_names

  # Create row names for each pairwise comparison
  contrast_names <- character(n_contrasts)

  # Reference level is the first level of each factor
  primary_ref <- primary_levels[1]
  secondary_ref <- secondary_levels[1]

  for (i in seq_along(primary_pairs)) {
    primary1 <- primary_pairs[[i]][1]
    primary2 <- primary_pairs[[i]][2]
    contrast_names[i] <- paste0(
      primary1,
      " - ",
      primary2,
      " (",
      secondary_level,
      ")"
    )

    # Calculate coefficients for primary1 in given secondary level
    primary1_coefs <- rep(0, n_coefs)
    names(primary1_coefs) <- coef_names

    # Intercept always contributes (represents primary_ref + secondary_ref)
    primary1_coefs["(Intercept)"] <- 1

    # Add primary main effect if not reference primary
    if (primary1 != primary_ref) {
      primary1_main <- paste0(primary_var, primary1)
      if (primary1_main %in% coef_names) {
        primary1_coefs[primary1_main] <- 1
      }
    }

    # Add secondary main effect if not reference secondary
    if (secondary_level != secondary_ref) {
      secondary_main <- paste0(secondary_var, secondary_level)
      if (secondary_main %in% coef_names) {
        primary1_coefs[secondary_main] <- 1
      }
    }

    # Add interaction effect if both are non-reference
    if (primary1 != primary_ref && secondary_level != secondary_ref) {
      interaction_term <- paste0(
        primary_var,
        primary1,
        ":",
        secondary_var,
        secondary_level
      )
      if (interaction_term %in% coef_names) {
        primary1_coefs[interaction_term] <- 1
      }
    }

    # Calculate coefficients for primary2 in given secondary level
    primary2_coefs <- rep(0, n_coefs)
    names(primary2_coefs) <- coef_names

    # Intercept always contributes
    primary2_coefs["(Intercept)"] <- 1

    # Add primary main effect if not reference primary
    if (primary2 != primary_ref) {
      primary2_main <- paste0(primary_var, primary2)
      if (primary2_main %in% coef_names) {
        primary2_coefs[primary2_main] <- 1
      }
    }

    # Add secondary main effect if not reference secondary
    if (secondary_level != secondary_ref) {
      secondary_main <- paste0(secondary_var, secondary_level)
      if (secondary_main %in% coef_names) {
        primary2_coefs[secondary_main] <- 1
      }
    }

    # Add interaction effect if both are non-reference
    if (primary2 != primary_ref && secondary_level != secondary_ref) {
      interaction_term <- paste0(
        primary_var,
        primary2,
        ":",
        secondary_var,
        secondary_level
      )
      if (interaction_term %in% coef_names) {
        primary2_coefs[interaction_term] <- 1
      }
    }

    # Contrast is primary1 - primary2
    contrast_matrix[i, ] <- primary1_coefs - primary2_coefs
  }

  rownames(contrast_matrix) <- contrast_names
  contrast_matrix
}

#' Create pairwise contrasts for primary covariate with no interaction
#'
#' This function creates contrast matrices for testing all pairwise
#' comparisons between levels of a primary covariate when there is no
#' significant interaction with the secondary covariate. The model structure
#' is assumed to be y ~ primary + secondary + (1 | subject).
#'
#' @param model A fitted model object with additive effects (no interaction)
#' @param primary_var Character string specifying the primary variable name
#' @return A contrast matrix with rows for each pairwise primary comparison
#'   and columns for model coefficients
#' @export
create_pairwise_primary_main_effect_contrasts <- function(
  model,
  primary_var
) {
  # Get model coefficients and check inputs
  coef_names <- names(lme4::fixef(model))

  # Get primary levels from the model data
  primary_levels <- levels(model@frame[[primary_var]])

  # Create all pairwise combinations of primary levels
  primary_pairs <- utils::combn(primary_levels, 2, simplify = FALSE)
  n_contrasts <- length(primary_pairs)
  n_coefs <- length(coef_names)

  # Initialize contrast matrix
  contrast_matrix <- matrix(
    0,
    nrow = n_contrasts,
    ncol = n_coefs
  )
  colnames(contrast_matrix) <- coef_names

  # Create row names for each pairwise comparison
  contrast_names <- character(n_contrasts)

  # Reference level is the first level of the primary factor
  primary_ref <- primary_levels[1]

  for (i in seq_along(primary_pairs)) {
    primary1 <- primary_pairs[[i]][1]
    primary2 <- primary_pairs[[i]][2]
    contrast_names[i] <- paste0(primary1, " - ", primary2)

    # For main effects only, we just need the difference between
    # the primary covariate coefficients

    # Coefficient for primary1
    if (primary1 == primary_ref) {
      # Reference level contributes 0 to the contrast
      primary1_coef <- 0
    } else {
      # Non-reference level contributes its main effect coefficient
      primary1_main <- paste0(primary_var, primary1)
      if (primary1_main %in% coef_names) {
        contrast_matrix[i, primary1_main] <- 1
      } else {
        warning(paste(
          "Main effect term",
          primary1_main,
          "not found in model coefficients"
        ))
      }
    }

    # Coefficient for primary2
    if (primary2 == primary_ref) {
      # Reference level contributes 0 to the contrast (already initialized)
    } else {
      # Non-reference level contributes negative of its main effect coefficient
      primary2_main <- paste0(primary_var, primary2)
      if (primary2_main %in% coef_names) {
        contrast_matrix[i, primary2_main] <- -1
      } else {
        warning(paste(
          "Main effect term",
          primary2_main,
          "not found in model coefficients"
        ))
      }
    }
  }

  rownames(contrast_matrix) <- contrast_names
  contrast_matrix
}
