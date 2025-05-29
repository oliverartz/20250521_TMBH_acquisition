## -----------------------------------------------------------------------------
## Purpose of script: Define parameters and functions for plotting
##
## Author: Oliver Artz
## Date Created: May 26, 2025
## Email: artzo@mskcc.org
##
## -----------------------------------------------------------------------------

# define general theme for plotting
theme_jco <- function() {
  theme_minimal(base_size = 10) +
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
      #panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      
      # customize legend
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.background = element_blank(),
      
      # customize facets
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey95", colour = "black", size = 0.5),
      
      # title adjustments
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      
      # spacious margins
      plot.margin = margin(10, 10, 10, 10)
    )
}

# function for plotting of variant types
plot_variant_type <- function(type, acquired_status) {
  
  # DEBUG ----
  if (FALSE){
    type <- "Indel"
    acquired_status <- "N"
  }
  # END DEBUG ----
  
  # get type to plot
  type_of_interest <- type
  
  # filter df
  df_plot <- merge_subtype %>%
    filter(type == type_of_interest) %>% 
    filter(acquired_tmbh == acquired_status)
  
  # make title
  n_total_title <- df_plot$n %>% sum()
  frac_total <- n_total_title / unique(df_plot$n_total)
  title_text <- paste0(type_of_interest, " (n=", n_total_title,", ", round(frac_total*100, 0), "%)")
  
  # plot
  df_plot %>%
    ggplot(aes(x = type, y = frac_subtype, fill = subtype)) +
    geom_bar(stat = "identity", color = "black") +
    # geom_text(aes(label = round(frac_subtype, 2), y = frac_subtype),
    #   position = position_stack(vjust = 0.5)) + 
    theme_jco() +
    theme(
      legend.position = "right",
      panel.grid = element_blank()) +
    labs(
      title = title_text,
      x = "",
      y = paste0("Frequency of ", type, " type"),
      fill = paste0(type, " Type")) +
    scale_fill_jco() +
    facet_wrap(~ acquired_tmbh, 
               nrow = 2)
}
