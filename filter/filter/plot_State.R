plot_State <- function(x, output_file = NULL) {
    # Convert matrices to data frames for easier plotting
    x_df <- as.data.frame(x)
    
    # Add time column
    x_df$time <- seq(0, T, length.out = nrow(x_df))
    
    # Reshape data frames for plotting
    x_melt <- reshape2::melt(x_df, id.vars = "time", variable.name = "Dimension", value.name = "Value")
    
    # Define a color palette that can accommodate more than 6 dimensions
    colors <- colorRampPalette(c('blue', 'red', 'green', 'purple', 'orange', 'brown', 'pink', 'cyan'))(length(unique(x_melt$Dimension)))
    plots <- list()
    for (i in unique(x_melt$Dimension)) {
      color <- colors[(as.numeric(gsub('V', '', i)) - 1) %% length(colors) + 1]
      p <- ggplot(subset(x_melt, Dimension == i), aes(x = time, y = Value)) +
        geom_line(color = color) +
        labs(title = paste0("x_", i), x = "Time", y = "Value") +
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
