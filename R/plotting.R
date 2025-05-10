#' Plot X to Y causal effect across different R values
#'
#' @param df_cleaned Cleaned data frame from process_simulation_results
#' @param save_plot Whether to save the plot to a file (default: FALSE)
#' @param filename Filename to save plot (default: "Bi_TSLS_Strength_XY_lines.png")
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 6)
#'
#' @return ggplot object
#' @export
#'
#' @import ggplot2
plot_bi_tsls_strength_XY <- function(df_cleaned, save_plot = FALSE,
                                     filename = "Bi_TSLS_Strength_XY_lines.png",
                                     width = 10, height = 6) {
  df_cleaned$R_w <- as.factor(df_cleaned$R_w)
  title_XY <- bquote(hat(beta)[X %->% Y] ~ "across different values of" ~ R[w] ~ "and" ~ -R[z])

  # Create base theme
  base_theme <- theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "#F5F5F5"),
      panel.grid.major = element_line(color = "lightgray"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "none",
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16)
    )

  # Calculate summary statistics for each group
  summary_data <- df_cleaned %>%
    dplyr::group_by(R_w) %>%
    dplyr::summarize(
      mean_XY = mean(beta_XY, na.rm = TRUE),
      min_XY = min(beta_XY, na.rm = TRUE),
      max_XY = max(beta_XY, na.rm = TRUE)
    )

  # Plot for X → Y with vertical lines and points
  p_XY <- ggplot(summary_data, aes(x = R_w)) +
    geom_errorbar(aes(ymin = min_XY, ymax = max_XY),
                  width = 0.2, color = "black", size = 0.7) +
    geom_point(aes(y = mean_XY), color = "#E74C3C", size = 3) +
    scale_x_discrete(
      name = bquote("Sensitivity analysis for different values of" ~ R[w] ~ "and" ~ -R[z]),
      labels = as.character(as.numeric(levels(summary_data$R_w)))
    ) +
    scale_y_continuous(
      name = "Point and interval estimates",
      labels = scales::label_number(accuracy = 0.01)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    base_theme +
    labs(title = title_XY)

  # Save plot if requested
  if (save_plot) {
    ggsave(filename, plot = p_XY, width = width, height = height,
           units = "in", dpi = 300)
  }

  return(p_XY)
}

#' Plot Y to X causal effect across different R values
#'
#' @param df_cleaned Cleaned data frame from process_simulation_results
#' @param save_plot Whether to save the plot to a file (default: FALSE)
#' @param filename Filename to save plot (default: "Bi_TSLS_Strength_YX_lines.png")
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 6)
#'
#' @return ggplot object
#' @export
#'
#' @import ggplot2
plot_bi_tsls_strength_YX <- function(df_cleaned, save_plot = FALSE,
                                     filename = "Bi_TSLS_Strength_YX_lines.png",
                                     width = 10, height = 6) {
  df_cleaned$R_w <- as.factor(df_cleaned$R_w)
  title_YX <- bquote(hat(beta)[Y %->% X] ~ "across different values of" ~ R[w] ~ "and" ~ -R[z])

  # Create base theme
  base_theme <- theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "#F5F5F5"),
      panel.grid.major = element_line(color = "lightgray"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "none",
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16)
    )

  # Calculate summary statistics for each group
  summary_data <- df_cleaned %>%
    dplyr::group_by(R_w) %>%
    dplyr::summarize(
      mean_YX = mean(beta_YX, na.rm = TRUE),
      min_YX = min(beta_YX, na.rm = TRUE),
      max_YX = max(beta_YX, na.rm = TRUE)
    )

  # Plot for Y → X with vertical lines and points
  p_YX <- ggplot(summary_data, aes(x = R_w)) +
    geom_errorbar(aes(ymin = min_YX, ymax = max_YX),
                  width = 0.2, color = "black", size = 0.7) +
    geom_point(aes(y = mean_YX), color = "#E74C3C", size = 3) +
    scale_x_discrete(
      name = bquote("Sensitivity analysis for different values of" ~ R[w] ~ "and" ~ -R[z]),
      labels = as.character(as.numeric(levels(summary_data$R_w)))
    ) +
    scale_y_continuous(
      name = "Point and interval estimates",
      labels = scales::label_number(accuracy = 0.01)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    base_theme +
    labs(title = title_YX)

  # Save plot if requested
  if (save_plot) {
    ggsave(filename, plot = p_YX, width = width, height = height,
           units = "in", dpi = 300)
  }

  return(p_YX)
}

#' Create combined plots of XY and YX effects
#'
#' @param df_cleaned Cleaned data frame with simulation results
#' @param save_plot Whether to save the plots (default: FALSE)
#' @param filename_prefix Prefix for filenames if saving (default: "Bi_TSLS_")
#'
#' @return List with XY and YX plots
#' @export
create_bi_tsls_plots <- function(df_cleaned, save_plot = FALSE, filename_prefix = "Bi_TSLS_") {
  # Create the plots
  xy_plot <- plot_bi_tsls_strength_XY(
    df_cleaned,
    save_plot = save_plot,
    filename = paste0(filename_prefix, "Strength_XY_lines.png")
  )

  yx_plot <- plot_bi_tsls_strength_YX(
    df_cleaned,
    save_plot = save_plot,
    filename = paste0(filename_prefix, "Strength_YX_lines.png")
  )

  # Return both plots
  return(list(
    xy_plot = xy_plot,
    yx_plot = yx_plot
  ))
}
