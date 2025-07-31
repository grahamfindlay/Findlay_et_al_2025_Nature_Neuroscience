#' Color and shape utilities
#'
#' Functions and constants for consistent plotting colors and shapes.

#' Convert RGB values (0-1 range) to hex colors
#'
#' @param r Red value (0-1)
#' @param g Green value (0-1)
#' @param b Blue value (0-1)
#' @return Hex color string
#' @importFrom grDevices rgb
#' @export
rgb_to_hex <- function(r, g, b) {
  rgb(r, g, b, maxColorValue = 1)
}

#' Custom colors for experiments
#'
#' Named vector of colors for different experimental conditions.
#'
#' @export
experiment_colors <- c(
  "Novelty" = rgb_to_hex(0.4, 0.4, 0.4),
  "Locomotion" = rgb_to_hex(
    0.6509803921568628,
    0.4627450980392157,
    0.11372549019607843
  ),
  "Dual" = rgb_to_hex(
    0.9019607843137255,
    0.6705882352941176,
    0.00784313725490196
  )
)

#' Custom shapes for experiments
#'
#' Named vector of point shapes for different experimental conditions.
#'
#' @export
experiment_shapes <- c("Novelty" = 16, "Locomotion" = 17, "Dual" = 18)
