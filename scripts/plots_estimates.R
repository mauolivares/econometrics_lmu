#' @title Plots
#' 
#' @description A collection of functions to plot data from Monte Carlo exercise as a function of sample sizes to illustrate LLNs and CLTs.
#' 
#' @param sample_sizes Numeric. Different sample sizes.
#' @param estimates Matrix. Matrix containing the estimates from the Monte Carlo exercise.
#' 
#' @return Plots
#' 
#' @author Mauricio Olivares
#' @import ggplot2 dplyr


my_plots <- function(sample_sizes, estimates){
  # Convert results to a data frame
  results_df <- do.call(rbind, lapply(seq_along(sample_sizes), function(i) {
    data.frame(
      sample_size = sample_sizes[i],
      beta2 = estimates[i,]
    )
  }))
  
  # Plot the distribution of the estimated slopes
   density_plot <- ggplot(results_df, aes(x = beta2, fill = as.factor(sample_size))) +
    geom_density(alpha = 0.5) +
    facet_wrap(~sample_size, scales = "free") +
    labs(title = "Distribution of Slope Estimates",
         x = "Estimated Slope",
         y = "Density",
         fill = "Sample Size") +
    theme_minimal()
  
  # Plot the mean and variance of the estimated slopes
   summary_df <- results_df %>%
    group_by(sample_size) %>%
    summarise(
      mean_slope = mean(beta2),
      var_slope = var(beta2)
    )
    return(density_plot)
}

plot_asympt <- function(sample_sizes, estimates) {
  # Convert results to a data frame
  results_df <- do.call(rbind, lapply(seq_along(sample_sizes), function(i) {
    data.frame(
      sample_size = sample_sizes[i],
      beta2 = estimates[i, ]
    )
  }))
  
  # Plot the distribution of the estimated slopes
  density_plot <- ggplot(results_df, aes(x = beta2, fill = as.factor(sample_size))) +
    geom_density(alpha = 0.5) + # Adjust alpha for transparency
    labs(title = "Distribution of Slope Estimates",
         x = "Estimate",
         y = "Density",
         fill = "Sample Size") +
    theme_minimal() +
    theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  panel.grid.major = element_line(color = "grey90"),
  panel.grid.minor = element_line(color = "grey95"),
  axis.text = element_text(color = "black"),
  axis.title = element_text(color = "black"),
  plot.title = element_text(color = "black", face = "bold"),
  plot.subtitle = element_text(color = "black"),
  legend.background = element_rect(fill = "white", color = NA),
  legend.key = element_rect(fill = "white", color = NA),
  legend.text = element_text(color = "black"),
  legend.title = element_text(color = "black")
) 

  # Return the density plot
  return(density_plot)
}