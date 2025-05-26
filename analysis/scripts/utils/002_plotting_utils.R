## -----------------------------------------------------------------------------
## Purpose of script: Define parameters and functions for plotting
##
## Author: Oliver Artz
## Date Created: May 26, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

theme_jco <- function() {
  theme_minimal(base_size = 10, base_family = "Arial") +
    theme(
      # add light axis lines and ticks
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      
      # adjust axis text and title
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      
      # remove gridlines for a clean look
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      
      # customize legend
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.background = element_blank(),
      
      # customize facets
      strip.text = element_text(face = "bold",
                                size = 10),
      strip.background = element_rect(fill = "grey95", colour = "black", size = 0.5),
      
      # title adjustments
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      
      # spacious margins
      plot.margin = margin(10, 10, 10, 10)
    )
}
