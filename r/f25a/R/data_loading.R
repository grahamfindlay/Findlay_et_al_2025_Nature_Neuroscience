#' Data loading functions
#'
#' Functions to load various datasets from the package extdata directory.

#' Load measures data
#'
#' @param rev Logical, whether to load revised multicortical version
#' @return Data frame with measures data
#' @export
load_measures <- function(rev = FALSE) {
  if (rev) {
    path <- system.file(
      "extdata",
      "multicortical_condition_measures.pqt",
      package = "f25a"
    )
  } else {
    path <- system.file(
      "extdata",
      "condition_measures.pqt",
      package = "f25a"
    )
  }
  d <- arrow::read_parquet(path)
  exp_levels <- c("Novelty", "Locomotion", "Dual")
  d$experiment <- factor(d$experiment, levels = exp_levels)
  d$subject <- factor(d$subject)
  d
}

#' Load condition differences data
#'
#' @param rev Logical, whether to load revised multicortical version
#' @return Data frame with condition differences data
#' @export
load_condition_differences <- function(rev = FALSE) {
  if (rev) {
    path <- system.file(
      "extdata",
      "multicortical_condition_contrasts.pqt",
      package = "f25a"
    )
  } else {
    path <- system.file("extdata", "condition_contrasts.pqt", package = "f25a")
  }
  d <- arrow::read_parquet(path)
  exp_levels <- c("Novelty", "Locomotion", "Dual")
  d$experiment <- factor(d$experiment, levels = exp_levels)
  d$subject <- factor(d$subject)
  if (rev) {
    d$acronym <- factor(d$acronym)
    d$region <- factor(d$region)
  }
  d
}

#' Load measures for homeostasis testing
#'
#' @param rev Logical, whether to load revised multicortical version
#' @return Data frame filtered for homeostasis conditions
#' @export
load_measures_for_homeostasis_testing <- function(rev = FALSE) {
  d <- load_measures(rev = rev)

  condition_levels <- c(
    "Early.BSL.NREM",
    "Late.BSL.NREM",
    "Early.EXT.Wake",
    "Late.EXT.Wake",
    "Early.REC.NREM",
    "Late.REC.NREM"
  )
  d <- subset(d, d$condition %in% condition_levels)
  d$condition <- factor(d$condition, levels = condition_levels)

  d
}

#' Load SG/FG data by epoch type
#'
#' @return Data frame with SG/FG data
#' @export
load_sgfg_by_epoch_type <- function() {
  path <- system.file(
    "extdata",
    "AeryJones_sgfg_medians_by_epoch_type.pqt",
    package = "f25a"
  )
  d <- arrow::read_parquet(path)

  # Ensure proper level ordering, for contrast matrix coding
  exp_levels <- c("Novelty", "Locomotion", "Dual")
  d$Experiment <- factor(d$Experiment, levels = exp_levels)

  epoch_levels <- c("SPW.Wake", "No.SPW.Wake", "SPW.NREM", "No.SPW.NREM")
  d$`Epoch Type` <- factor(d[["Epoch.Type"]], levels = epoch_levels)

  d$Subject <- factor(d$Subject)

  d <- d[d$ROI == "CA1-slm", ]
  d
}

#' Load measures for extended wake change testing
#'
#' @return Data frame filtered for extended wake conditions
#' @export
load_measures_for_ewk_change_testing <- function() {
  d <- load_measures()

  condition_levels <- c("Early.EXT.Wake", "Late.EXT.Wake")
  d <- subset(d, d$condition %in% condition_levels)
  d$condition <- factor(d$condition, levels = condition_levels)

  d
}

#' Load measures for first novelty testing
#'
#' @param only_dual_subjects Logical, whether to include only dual experiment subjects
#' @return Data frame filtered for first novelty conditions
#' @export
load_measures_for_first_novelty_testing <- function(only_dual_subjects) {
  d <- load_measures()
  d <- droplevels(subset(d, d$condition %in% c("Early.NOD.Wake")))
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
  d
}

#' Load ripple frequency data
#'
#' @return Data frame with ripple frequency data
#' @export
load_ripple_frequency_data <- function() {
  path <- system.file(
    "extdata",
    "ripple_frequency_by_state.pqt",
    package = "f25a"
  )
  d <- arrow::read_parquet(path)

  exp_levels <- c("Novelty", "Locomotion", "Dual")
  d$experiment <- factor(d$experiment, levels = exp_levels)

  state_levels <- c("Wake", "NREM")
  d$state <- factor(d$state, state_levels)

  d$subject <- factor(d$subject)

  d
}

#' Load occupancy data
#'
#' @return Data frame with sleep period fractional occupancy data
#' @export
load_occupancy_data <- function() {
  # Eugene has already been excluded.
  path <- system.file(
    "extdata",
    "sleep_period_fractional_occupancy.pqt",
    package = "f25a"
  )
  d <- arrow::read_parquet(path)

  exp_levels <- c("Novelty", "Locomotion", "Dual")
  d$experiment <- factor(d$experiment, levels = exp_levels)

  day_levels <- c("Baseline", "Recovery")
  d$day <- factor(d$day, levels = day_levels)

  state_levels <- c("NREM", "IS", "REM", "Wake", "MA")
  d$state <- factor(d$state, levels = state_levels)

  d$subject <- factor(d$subject)

  d <- droplevels(subset(d, d$experiment %in% c("Novelty", "Locomotion")))
  d
}

#' Load recovery REM vs extended wake TDR data
#'
#' @return Data frame with recovery REM vs extended wake TDR data
#' @export
load_rec_rem_vs_ewk_tdr_data <- function() {
  path <- system.file("extdata", "rec_rem_vs_ewk_tdr.pqt", package = "f25a")
  dat <- arrow::read_parquet(path)

  dat <- dat[stats::complete.cases(dat), ] # Drop (intentionally) missing values

  exp_levels <- c("Novelty", "Locomotion", "Dual")
  dat$Experiment <- factor(dat$Experiment, levels = exp_levels)
  dat$Subject <- factor(dat$Subject)

  # Rename columns to use dot notation
  names(dat)[names(dat) == "REM fraction"] <- "REM.Fraction"
  names(dat)[
    names(dat) == "Total Hippocampal Theta:Delta"
  ] <- "Total.Hippocampal.ThetaDeltaRatio"

  dat
}

#' Load firing rates data
#'
#' @param response_var Character string specifying the response variable
#'   column name
#' @param region Character string specifying the brain region to filter for
#' @param cell_type Character string specifying the cell type to filter for
#' @param exclude_subjects_without_offs Logical, whether to exclude subjects
#'   without OFF firing rates (default: FALSE)
#' @return Data frame with firing rates data
#' @export
load_firing_rates <- function(
  response_var,
  region,
  cell_type,
  exclude_subjects_without_offs = FALSE
) {
  path <- system.file(
    "extdata",
    "firing_rates.pqt",
    package = "f25a"
  )
  d <- arrow::read_parquet(path)

  # Filter data
  d <- d[d$region == region, ]
  if (cell_type == "NA") {
    d <- d[is.na(d$petersen_cell_type), ]
  } else {
    d <- d[d$petersen_cell_type == cell_type, ]
  }
  if (exclude_subjects_without_offs) {
    d <- d[!is.na(d$ON.Firing.Rate), ]
  }
  d <- d[!is.na(d[[response_var]]), ]

  # Keep only the required columns
  d <- d[, c("subject", "cluster_id", "condition", response_var)]

  # Set factors
  condition_levels <- c(
    "Early.BSL.NREM",
    "Late.BSL.NREM",
    "Early.EXT.Wake",
    "Late.EXT.Wake",
    "Early.REC.NREM",
    "Late.REC.NREM"
  )
  d$condition <- factor(d$condition, levels = condition_levels)
  d$subject <- factor(d$subject)
  d
}
