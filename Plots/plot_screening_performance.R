# Data source: Tables 1-3 (rho = 0.2, Case 1, M1 imputation for baseline methods)
library(ggplot2)
library(dplyr)

plot_data <- data.frame(
  Model = c(rep("Linear", 10), rep("Poisson", 10), rep("Cox", 10)),
  Setting = rep(c(rep("Balanced", 5), rep("Unbalanced", 5)), 3),
  Method = rep(c("Log-Rank (Proposed)", "ADD-SIS", "DC-SIS", "KF", "MV-SIS"), 6),
  Median_MMS = c(
    17.0, 25.5, 27.0, 50.5, 28.0,  
    26.0, 52.0, 39.5, 58.5, 38.5,  

    10.0, 16.5, 44.0, 80.0, 53.0,  
    12.0, 18.0, 33.5, 59.2, 42.5,  

    18.0, 42.5, 39.0, 75.5, 50.5,  
    45.0, 80.0, 82.5, 112.5, 77.5  
  ),
  RSD = c(

    44.8, 114.6, 87.5, 133.0, 92.9,   
    89.7, 212.3, 141.2, 176.9, 134.3, 

    10.4, 32.1, 86.9, 120.0, 97.9,    
    14.2, 30.4, 78.7, 105.0, 82.5,    

    27.8, 81.7, 96.8, 125.4, 113.1,   
    56.7, 164.0, 132.3, 188.1, 138.1  
  )
)

y_axis_order <- rev(c("Log-Rank (Proposed)", "ADD-SIS", "DC-SIS", "MV-SIS", "KF"))
plot_data$Method <- factor(plot_data$Method, levels = y_axis_order)

legend_display_order <- c("Log-Rank (Proposed)", "ADD-SIS", "DC-SIS", "MV-SIS", "KF")

plot_data$Model <- factor(plot_data$Model, levels = c("Linear", "Poisson", "Cox"))

plot_data$Setting <- factor(plot_data$Setting, levels = c("Balanced", "Unbalanced"))


p <- ggplot(plot_data, aes(x = Median_MMS, y = Method)) +
  
  geom_segment(aes(xend = Median_MMS + RSD, yend = Method, color = Method), 
               linewidth = 1.2, alpha = 0.75) + 
  
  geom_point(aes(color = Method, shape = Method), size = 4.5, fill = "white", stroke = 1.5) +
  
  facet_grid(Setting ~ Model, scales = "free_x") +
  
  scale_color_manual(
    values = c(
      "Log-Rank (Proposed)" = "#D50000", 
      "ADD-SIS"             = "#1565C0", 
      "DC-SIS"              = "#2E7D32", 
      "MV-SIS"              = "#8E24AA", 
      "KF"                  = "#F57C00"  
    ),
    breaks = legend_display_order 
  ) +
  
  scale_shape_manual(
    values = c(19, 15, 17, 18, 1), 
    breaks = legend_display_order
  ) +
  
  scale_y_discrete(
    expand = expansion(mult = 0.1), 
    labels = c(
      "Log-Rank (Proposed)" = "Log-Rank\n(Proposed)",
      "ADD-SIS" = "ADD-SIS",
      "DC-SIS"  = "DC-SIS",
      "MV-SIS"  = "MV-SIS",
      "KF"      = "KF"
    )
  ) +
  
  theme_bw(base_size = 18, base_family = "serif") +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)), 
    axis.text.y  = element_text(size = 16, color = "black", face = "bold", lineheight = 0.9), 
    axis.text.x  = element_text(size = 15, color = "black", face = "bold"), 
    
    strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.8),
    
    strip.text.x = element_text(
      face = "bold", 
      size = 18,  
      margin = margin(t = 10, b = 10) 
    ),
    
    strip.text.y = element_text(
      face = "bold", 
      size = 18,  
      angle = 270,  
      margin = margin(t = 12, r = 12, b = 12, l = 12) 
    ),
    
    panel.spacing = unit(0.2, "lines"), 
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(), 
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.5, linetype = "dashed"),
    
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"), 
    legend.margin = margin(t = 5, b = 0),     
    legend.box.margin = margin(t = -2, b = 0), 
    legend.key.size = unit(1.2, "cm") 
  ) +
  
  labs(
    x = "Median Minimum Model Size (MMS) with Robust Standard Deviation (RSD)",
    y = NULL 
  )

print(p)

ggsave(filename = "Screening_Performance_2x3_LargeFont.pdf", 
       plot = p, 
       width = 13,   
       height = 8.5,     
       device = "pdf")