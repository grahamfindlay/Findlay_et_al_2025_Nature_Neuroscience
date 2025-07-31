#' Contrast matrices for sleep homeostasis and extended wake testing

#' Tests for homeostatic effects in novelty experiments, in the presence of an
#' interction between experiment and condition.
#'
#' @format Matrix (5x18) with contrasts:
#'     E.REC-L.BSL, L.REC-E.REC, L.EWK-E.EWK, E.BSL-L.BSL, E.REC-E.BSL
#' @export
homeostasis_novelty_matrix <- rbind(
  "E.REC - L.BSL" = c(0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "L.REC - E.REC" = c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "L.EWK - E.EWK" = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "E.BSL - L.BSL" = c(0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "E.REC - E.BSL" = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)

#' Tests for homeostatic effects in locomotion experiments, in the presence of
#' an interction between experiment and condition.
#'
#' @format Matrix (5x18) with contrasts:
#'     E.REC-L.BSL, L.REC-E.REC, L.EWK-E.EWK, E.BSL-L.BSL, E.REC-E.BSL
#' @export
homeostasis_locomotion_matrix <- rbind(
  "E.REC - L.BSL" = c(0, -1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0),
  "L.REC - E.REC" = c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0),
  "L.EWK - E.EWK" = c(0, 0, -1, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0),
  "E.BSL - L.BSL" = c(0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "E.REC - E.BSL" = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
)

#' Tests for homeostatic effects in dual experiments, in the presence of an
#' interction between experiment and condition.
#'
#' @format Matrix (5x18) with contrasts:
#'     E.REC-L.BSL, L.REC-E.REC, L.EWK-E.EWK, E.BSL-L.BSL, E.REC-E.BSL
#' @export
homeostasis_dual_matrix <- rbind(
  "E.REC - L.BSL" = c(0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0),
  "L.REC - E.REC" = c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1),
  "L.EWK - E.EWK" = c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0),
  "E.BSL - L.BSL" = c(0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0),
  "E.REC - E.BSL" = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
)

#' Tests for homeostatic effects when there is no interaction between experiment
#' and condition, and a significant main effect of condition.
#'
#' @format Matrix (5x8) with contrasts:
#'     E.REC-L.BSL, L.REC-E.REC, L.EWK-E.EWK, E.BSL-L.BSL, E.REC-E.BSL
#' @export
homeostasis_main_effect_matrix <- rbind(
  "E.REC - L.BSL" = c(0, -1, 0, 0, 1, 0, 0, 0),
  "L.REC - E.REC" = c(0, 0, 0, 0, -1, 1, 0, 0),
  "L.EWK - E.EWK" = c(0, 0, -1, 1, 0, 0, 0, 0),
  "E.BSL - L.BSL" = c(0, -1, 0, 0, 0, 0, 0, 0),
  "E.REC - E.BSL" = c(0, 0, 0, 0, 1, 0, 0, 0)
)

# Extended Wake change matrices -----------------------------------------------

#' Tests for extended wake change in novelty experiments, in the presence of an
#' interction between experiment and condition.
#'
#' @format Matrix (1x6) with contrast: L.EWK-E.EWK
#' @export
ewk_change_novelty_matrix <-
  rbind("L.EWK - E.EWK" = c(0, 1, 0, 0, 0, 0))

#' Tests for extended wake change in locomotion experiments, in the presence of
#' an interction between experiment and condition.
#'
#' @format Matrix (1x6) with contrast: L.EWK-E.EWK
#' @export
ewk_change_locomotion_matrix <-
  rbind("L.EWK - E.EWK" = c(0, 1, 0, 0, 1, 0))

#' Tests for extended wake change in dual experiments, in the presence of an
#' interction between experiment and condition.
#'
#' @format Matrix (1x6) with contrast: L.EWK-E.EWK
#' @export
ewk_change_dual_matrix <-
  rbind("L.EWK - E.EWK" = c(0, 1, 0, 0, 0, 1))

#' Tests for extended wake change when there is no interaction between
#' experiment and condition, and a significant main effect of condition.
#'
#' @format Matrix (1x4) with contrast: L.EWK-E.EWK
#' @export
ewk_change_main_effect_matrix <-
  rbind("L.EWK - E.EWK" = c(0, 1, 0, 0))

# Correlation matrices --------------------------------------------------------

#' Test for significant correlation between a continuous response variable and
#' a continuous covariate, in the presence of an interaction between experiment
#' and covariate.
#'
#' @format Matrix (3x6) with experiment contrasts: Novelty, Locomotion, Dual
#' @export
correlation_interaction_matrix <- rbind(
  "Novelty" = c(0, 1, 0, 0, 0, 0),
  "Locomotion" = c(0, 1, 0, 0, 1, 0),
  "Dual" = c(0, 1, 0, 0, 0, 1)
)

# Occupancy matrices ----------------------------------------------------------

#' Tests for significant differences in occupancy between COW and NOD, by
#' sleep state.
#'
#' @format Matrix (5x10) with COW-NOD contrasts: NREM, IS, REM, Wake, MA
#' @export
day_occupancy_contrast_matrix <- rbind(
  "NREM: COW - NOD" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  "IS: COW - NOD" = c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
  "REM: COW - NOD" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0),
  "Wake: COW - NOD" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0),
  "MA: COW - NOD" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1)
)

#' Tests for significant differences in occupancy between REC and BSL, by
#' sleep state.
#'
#' @format Matrix (5x10) with REC-BSL contrasts: NREM, IS, REM, Wake, MA
#' @export
exp_occupancy_contrast_matrix <- rbind(
  "NREM: REC - BSL" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  "IS: REC - BSL" = c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
  "REM: REC - BSL" = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0),
  "Wake: REC - BSL" = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0),
  "MA: REC - BSL" = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1)
)
