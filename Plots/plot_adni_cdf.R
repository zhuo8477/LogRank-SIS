# Prerequisites:
# Run the ADNI data analysis script first to load:
#   l, r                  - interval endpoints (months to dementia)
#   X_mat                 - SNP binary predictor matrix
#   snp_names             - SNP identifiers
#   rank_lr_adni          - Log-rank SIS ranks (genome-wide)
#   rank_add_adni         - ADD-SIS ranks
#   turnbull_em_npmle     - NPMLE EM algorithm function
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggtext)
library(patchwork)


pal_lines <- c(
  "Ref allele (X = 0)" = "#2C7BB6", 
  "Alt allele (X = 1)" = "#D7191C"
)
pal_types <- c(
  "Ref allele (X = 0)" = "solid", 
  "Alt allele (X = 1)" = "dashed"
)

get_npmle_cdf <- function(l_sub, r_sub, supp_all) {
  p <- tryCatch(
    turnbull_em_npmle(l_sub, r_sub, supp_all),
    error = function(e) rep(0, length(supp_all))
  )
  cumsum(p)
}

supp_global <- sort(unique(c(l[is.finite(l) & l > 0], r[is.finite(r)])))
x_max_global <- max(supp_global, na.rm = TRUE)

base_theme <- theme_classic(base_size = 15) +
  theme(
    axis.line.x         = element_line(colour = "black", linewidth = 0.8),
    axis.line.y         = element_line(colour = "black", linewidth = 0.8),
    axis.ticks          = element_line(colour = "black", linewidth = 0.8),
    axis.text           = element_text(size = 14, colour = "black"),
    axis.title          = element_text(size = 16, face = "bold", colour = "black"),
    
    panel.grid.major.y  = element_line(colour = "#e5e5e5", linewidth = 0.5, linetype = "dashed"),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.background    = element_rect(fill = "white", colour = NA),
    plot.background     = element_rect(fill = "white", colour = NA),
    
    plot.title          = element_markdown(size = 16, face = "bold", margin = margin(b = 12)),
    
    legend.key.width    = unit(4.5, "cm"),
    legend.key.height   = unit(0.6, "cm"),
    legend.background   = element_blank(),
    legend.box.background = element_blank()
  )


make_npmle_panel <- function(snp_name, panel_lab, main_title, rank_lr, rank_add) {
  
  x_k  <- X_mat[, snp_name]
  idx0 <- which(x_k == 0)
  idx1 <- which(x_k == 1)
  n0   <- length(idx0)
  n1   <- length(idx1)
  
  F0 <- get_npmle_cdf(l[idx0], r[idx0], supp_global)
  F1 <- get_npmle_cdf(l[idx1], r[idx1], supp_global)
  
  t_vec  <- c(0, supp_global, x_max_global)
  F0_vec <- c(0, F0, tail(F0, 1))
  F1_vec <- c(0, F1, tail(F1, 1))
  
  df_plot <- bind_rows(
    tibble(Time = t_vec, CDF = F0_vec, Group = "Ref allele (X = 0)"),
    tibble(Time = t_vec, CDF = F1_vec, Group = "Alt allele (X = 1)")
  ) %>%
    mutate(Group = factor(Group, levels = c("Ref allele (X = 0)", "Alt allele (X = 1)")))
  
  maf_v <- round(min(mean(x_k), 1 - mean(x_k)), 4)
  
  title_md <- sprintf("%s %s", panel_lab, main_title)
  
  anno_text <- sprintf(
    "<span style='color:black;'>**SNP: %s**</span><br>
     <span style='color:#333333; font-size:12pt;'>MAF = %.4f</span><br>
     <span style='color:#333333; font-size:12pt;'>*n* = %d vs %d</span><br>
     <span style='color:#333333; font-size:12pt;'>LR rank: %d &nbsp;&nbsp;|&nbsp;&nbsp; ADD rank: %d</span>",
    snp_name, maf_v, n0, n1, rank_lr, rank_add
  )
  
  p <- ggplot(df_plot, aes(x = Time, y = CDF, color = Group, linetype = Group)) +
    geom_step(linewidth = 1.3, direction = "hv") + 
    scale_color_manual(values = pal_lines) +
    scale_linetype_manual(values = pal_types) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03)), breaks = pretty_breaks(n = 6)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title = title_md,
      x = "Time to dementia (months)",
      y = expression(bolditalic(F)(bolditalic(t)))
    ) +
    
    annotate(
      "richtext", 
      x = x_max_global * 0.98, 
      y = 0.04, 
      label = anno_text, 
      hjust = 1, vjust = 0, 
      size = 4.8,                        
      fill = alpha("white", 0.85),       
      label.color = "#cccccc",           
      label.padding = unit(0.5, "lines"),
      label.r = unit(0.15, "lines")      
    ) +
    base_theme
  
  return(p)
}


get_rk <- function(snp) {
  idx <- which(snp_names == snp)
  c(LR = rank_lr_adni[idx], ADD = rank_add_adni[idx])
}

rk_a <- get_rk("rs12941385_A"); rk_b <- get_rk("rs4235426_G")
rk_c <- get_rk("rs1895140_G");  rk_d <- get_rk("rs1947078_A")

p_a <- make_npmle_panel("rs12941385_A", "(a)", "Consensus signal", rk_a["LR"], rk_a["ADD"])
p_b <- make_npmle_panel("rs4235426_G", "(b)", "Log-rank SIS unique", rk_b["LR"], rk_b["ADD"])
p_c <- make_npmle_panel("rs1895140_G", "(c)", "ADD-SIS unique (crossing hazards)", rk_c["LR"], rk_c["ADD"])
p_d <- make_npmle_panel("rs1947078_A", "(d)", "Null baseline", rk_d["LR"], rk_d["ADD"])

p_a <- p_a + theme(axis.title.x = element_blank())
p_b <- p_b + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p_d <- p_d + theme(axis.title.y = element_blank())

combined_plot <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(guides = "collect") & 
  theme(
    legend.position   = "bottom",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 15, face = "bold"),
    legend.margin     = margin(t = 5, b = 10),
    plot.margin       = margin(15, 15, 10, 15)
  ) &
  plot_annotation(
    title = "NPMLE cumulative incidence F(t) stratified by SNP genotype group (ADNI)",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 20))
    )
  )


ggsave(
  filename = "Fig_NPMLE_CDF_ggplot_Max.pdf", 
  plot     = combined_plot, 
  width    = 12.5,
  height   = 10.5, 
  device   = cairo_pdf, 
  dpi      = 600
)
