plot_Obser <- function(y_tau, output_file = NULL) {
    # Convert matrices to data frames for easier plotting
    y_tau_df <- as.data.frame(y_tau)
    
    # Add time column
    y_tau_df$time <- seq(0, T, length.out = nrow(y_tau_df))
    
    # Reshape data frames for plotting
    y_tau_melt <- reshape2::melt(y_tau_df, id.vars = "time", variable.name = "Dimension", value.name = "Value")
    
    # Define a color palette that can accommodate more than 6 dimensions
    colors <- colorRampPalette(c('blue', 'red', 'green', 'purple', 'orange', 'brown', 'pink', 'cyan'))(length(unique(y_tau_melt$Dimension)))
    plots <- list()
    for (i in unique(y_tau_melt$Dimension)) {
      color <- colors[(as.numeric(gsub('V', '', i)) - 1) %% length(colors) + 1]
      p <- ggplot(subset(y_tau_melt, Dimension == i), aes(x = time, y = Value)) +
        geom_line(color = color) +
        labs(title = paste0("y_", i), x = "Time", y = "Value") +
        theme_minimal() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = 'black', fill = NA, linewidth = 1))
      plots[[length(plots) + 1]] <- p
    }
    
    # Arrange the plots dynamically for better visualization
    n <- length(plots)
    if (n <= 3) {
      ncol <- n
      nrow <- 1
    } else {
      ncol <- ceiling(sqrt(n))
      nrow <- ceiling(n / ncol)
    }
    grid_plot <- do.call(grid.arrange, c(plots, ncol = ncol, nrow = nrow))
    
    # Save the plot if output_file is provided
    if (!is.null(output_file)) {
      ggsave(output_file, grid_plot, width = 12, height = ceiling(length(plots) / 3) * 4)
    }
  }
