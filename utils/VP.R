library(ggplot2)

# Create volcano plot
volcano_plot <- function(data) {
  # Calculate -log10(padj) for y-axis
  data$neg_log_padj <- -log10(data$padj)
  
  # Create color variable based on significance criteria
  data$color <- ifelse(abs(data$log2FoldChange) > 1 & data$padj < 0.05, "red", "black")
  
  # Create the plot
  p <- ggplot(data, aes(x = log2FoldChange, y = neg_log_padj)) +
    geom_point(aes(size = log2(baseMean), color = color), alpha = 0.7) +
    scale_color_identity() +  # Use the actual colors specified
    
    # Add threshold lines
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    
    # Labels and theme
    labs(
      title = "Volcano Plot",
      x = "Log2 Fold Change (logFC)",
      y = "-Log10(adjusted p-value)",
      size = "Base Mean"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.position = "right"
    ) +
    
    # Adjust point size range for better visualization
    scale_size_continuous(range = c(0.5, 3), guide = guide_legend(override.aes = list(color = "black")))
  
  return(p)
}
