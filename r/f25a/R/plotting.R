#' Plot a correlation with experiment interaction
#'
#' Creates a plot showing correlation results when there is a significant
#' interaction between a primary covariate and experiment.
#'
#' @param result Result object from correlation test with significant
#'   interaction
#' @param response_var Name of response variable (y-axis)
#' @param primary_covariate Name of primary covariate (x-axis)
#' @param fig_dir Directory to save figures to. If NULL, figures are not saved.
#' @return ggplot2 plot object
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon labs theme_minimal ggsave
#' @export
plot_correlation_with_experiment_interaction <- function(
  result,
  response_var,
  primary_covariate,
  fig_dir = NULL
) {
  # Extract covariate names from the model
  secondary_covariate <- "experiment"

  interaction_predictions <- modelbased::estimate_expectation(
    result$models$full,
    by = c(primary_covariate, secondary_covariate),
    length = 100
  )

  # Determine which experiments have significant or trending slopes
  posthoc_summary <- summary(result$interaction$posthoc$glht)
  posthoc_pvals <- posthoc_summary$test$pvalues

  # Include experiments with p < 0.1 (both significant and trending)
  relevant_experiments <- levels(result$data$experiment)[
    which(posthoc_pvals < 0.1)
  ]

  # Create a mapping of experiment to line type
  experiment_linetypes <- ifelse(
    posthoc_pvals < 0.05,
    "solid",
    "dashed"
  )
  names(experiment_linetypes) <- levels(result$data$experiment)

  # Filter predictions to include significant and trending experiments
  sig_predictions <- interaction_predictions %>%
    dplyr::filter(get(secondary_covariate) %in% relevant_experiments) %>%
    dplyr::mutate(
      linetype = experiment_linetypes[get(secondary_covariate)]
    )

  # Plot significant and trending slopes
  p <-
    ggplot2::ggplot(
      sig_predictions,
      ggplot2::aes(
        x = get(primary_covariate),
        y = Predicted,
        color = get(secondary_covariate),
        linetype = linetype
      )
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = CI_low,
        ymax = CI_high,
        fill = get(secondary_covariate)
      ),
      alpha = 0.2,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 1.2, alpha = 0.5) +
    ggplot2::geom_point(
      data = result$data,
      ggplot2::aes(
        x = get(primary_covariate),
        y = get(response_var),
        color = get(secondary_covariate),
        shape = get(secondary_covariate)
      ),
      size = 5,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(values = experiment_colors) +
    ggplot2::scale_fill_manual(values = experiment_colors) +
    ggplot2::scale_shape_manual(values = experiment_shapes) +
    ggplot2::scale_linetype_identity() +
    ggplot2::labs(
      x = primary_covariate,
      y = response_var,
      color = secondary_covariate,
      shape = secondary_covariate
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      aspect.ratio = 1
    )

  # Create version without axis labels or tick labels for saving
  # Note: Updated to use a configurable path instead of hard-coded here::here
  p_clean <- p +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  if (!is.null(fig_dir)) {
    plotname <- paste0(response_var, "_vs_", primary_covariate, ".png")
    plotpath <- file.path(fig_dir, plotname)
    ggplot2::ggsave(
      plotpath,
      p_clean,
      bg = "white",
      width = 5,
      height = 5,
      create.dir = TRUE
    )
  }

  p
}

#' Plot condition correlation without interaction
#'
#' Creates a plot showing correlation results when there is no significant
#' interaction with experiment.
#'
#' @param result Result object from correlation test without significant
#'   interaction
#' @param response_var Name of response variable (y-axis)
#' @param primary_covariate Name of primary covariate (x-axis)
#' @param fig_dir Directory to save figures to. If NULL, figures are not saved.
#' @return ggplot2 plot object
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon labs theme_minimal ggsave
#' @export
plot_correlation_without_experiment_interaction <- function(
  result,
  response_var,
  primary_covariate,
  fig_dir = NULL
) {
  # Extract covariate names from the model
  secondary_covariate <- "experiment"

  # Generate main effect predictions
  main_effect_predictions <- modelbased::estimate_expectation(
    result$models$reduced,
    by = primary_covariate,
    length = 100
  )

  # Determine line type based on main effect p-value
  show_main_effect <- result$main_effect$pval < 0.1
  main_effect_linetype <- if (result$main_effect$pval < 0.05) {
    "solid"
  } else {
    "dashed"
  }

  # Plot main effect as solid black line
  p <-
    ggplot2::ggplot(
      main_effect_predictions,
      ggplot2::aes(
        x = get(primary_covariate),
        y = Predicted
      )
    ) +
    {
      if (show_main_effect) {
        ggplot2::geom_ribbon(
          ggplot2::aes(
            ymin = CI_low,
            ymax = CI_high
          ),
          alpha = 0.2,
          fill = "black"
        )
      }
    } +
    {
      if (show_main_effect) {
        ggplot2::geom_line(
          color = "black",
          linewidth = 1.2,
          linetype = main_effect_linetype
        )
      }
    } +
    ggplot2::geom_point(
      data = result$data,
      ggplot2::aes(
        x = get(primary_covariate),
        y = get(response_var),
        color = get(secondary_covariate),
        shape = get(secondary_covariate)
      ),
      size = 5,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(values = experiment_colors) +
    ggplot2::scale_shape_manual(values = experiment_shapes) +
    ggplot2::labs(
      x = primary_covariate,
      y = response_var,
      color = secondary_covariate,
      shape = secondary_covariate
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      aspect.ratio = 1
    )

  # Create version without axis labels or tick labels for saving
  p_clean <- p +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  if (!is.null(fig_dir)) {
    plotname <- paste0(response_var, "_vs_", primary_covariate, ".png")
    plotpath <- file.path(fig_dir, plotname)
    ggplot2::ggsave(
      plotpath,
      p_clean,
      bg = "white",
      width = 5,
      height = 5,
      create.dir = TRUE
    )
  }

  p
}
