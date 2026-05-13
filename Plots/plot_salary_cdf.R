# Prerequisites:
# Run the job advertisement data analysis script first to load:
#   l_fig, r_fig, s_fig   - interval endpoints and binary predictors
#   word_names            - keyword labels
#   GLOBAL_MAX_X          - x-axis upper bound
#   words_cfg             - panel configuration list
#   pull_scores_and_ranks - helper function for ADD-SIS and Log-rank scores
#   emicm_to_df, EMICM    - NPMLE helper functions
library(ggplot2)
library(patchwork)

make_panel <- function(cfg) {
  word_str    <- cfg$word
  word_label  <- cfg$label
  tag         <- cfg$tag
  
  idx <- which(word_names == word_str)[1]
  
  if(is.na(idx)) {
    return(ggplot() + theme_void() + labs(title = paste(tag, word_label, "[Word Not Found]")))
  }
  
  x_k   <- as.integer(s_fig[[idx]] == TRUE | s_fig[[idx]] == "TRUE" | s_fig[[idx]] == 1)
  n_pos <- sum(x_k == 1)
  n_neg <- sum(x_k == 0)
  
  if (n_pos < 2 || n_neg < 2) {
    return(ggplot() + theme_void() + labs(title = paste(tag, word_label, "[Insufficient Data]")))
  }
  
  em_pos <- suppressMessages(suppressWarnings(EMICM(cbind(l_fig[x_k == 1], r_fig[x_k == 1]))))
  em_neg <- suppressMessages(suppressWarnings(EMICM(cbind(l_fig[x_k == 0], r_fig[x_k == 0]))))
  
  lbl_pos <- sprintf("With keyword")
  lbl_neg <- sprintf("Without keyword (Baseline)")
  
  df_all <- rbind(emicm_to_df(em_pos, lbl_pos), emicm_to_df(em_neg, lbl_neg))
  
  d_val <- if (exists("d")) d else 60
  info  <- pull_scores_and_ranks(idx)
  
  ann_txt <- if (!is.na(info$add_sc) && !is.na(info$lr_sc)) {
    sprintf("ADD-SIS : %9.3f (%s)\nLog-rank: %9.3f (%s)", 
            info$add_sc, fmt_rank(info$add_rk, d_val), 
            info$lr_sc,  fmt_rank(info$lr_rk, d_val))
  } else {
    "Score/Rank unavailable"
  }
  
  pal  <- setNames(c("#B2182B", "#2166AC"), c(lbl_pos, lbl_neg)) 
  ltys <- setNames(c("solid",   "dashed"),  c(lbl_pos, lbl_neg))
  szs  <- setNames(c(1.4,       1.1),       c(lbl_pos, lbl_neg)) 
  
  p <- ggplot(df_all, aes(x = x, y = cdf, colour = group, linetype = group, linewidth = group)) +
    geom_step(direction = "hv", alpha = 0.9) +
    annotate("text", x = GLOBAL_MAX_X * 0.96, y = 0.05, label = ann_txt,
             hjust = 1, vjust = 0, size = 6.5, family = "mono", fontface = "bold", 
             colour = "#222222", lineheight = 1.1) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_linetype_manual(values = ltys, name = NULL) +
    scale_linewidth_manual(values = szs, name = NULL) +
    scale_x_continuous(
      name   = "Posted Salary (yuan / month)",
      breaks = seq(0, 40000, by = 10000),
      labels = c("0", "10k", "20k", "30k", "40k"), expand = c(0, 0)
    ) +
    scale_y_continuous(
      name   = "Cumulative Probability",
      breaks = seq(0, 1, by = 0.25),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"), expand = c(0, 0)
    ) +
    coord_cartesian(xlim = c(0, GLOBAL_MAX_X), ylim = c(0, 1.05)) +
    labs(title = paste(tag, word_label)) +
    theme_bw(base_size = 20) +
    theme(
      plot.title        = element_text(size = 22, face = "bold", margin = margin(b = 5)),
      axis.title        = element_text(size = 20, face = "bold"),
      axis.text         = element_text(size = 18, face = "bold", colour = "black"),
      legend.position   = "bottom", 
      legend.key.width  = unit(3, "cm"),
      legend.text       = element_text(size = 18, face = "bold"),
      panel.grid.minor  = element_blank(), 
      panel.grid.major.x= element_blank(),
      panel.grid.major.y= element_line(colour = "#E8E8E8", linewidth = 0.4),
      panel.border      = element_rect(colour = "black", linewidth = 1.2), 
      plot.margin       = margin(5, 5, 5, 5)
    )
  
  return(p)
}

panels <- lapply(words_cfg, make_panel)
names(panels) <- sapply(words_cfg, `[[`, "tag")

fig_out <- (panels[["(a)"]] | panels[["(b)"]]) /
  (panels[["(c)"]] | panels[["(d)"]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box.margin = margin(t = 5, b = 5))

print(fig_out)

ggsave(
  filename = "salary_analysis_plot.pdf", 
  plot = fig_out, 
  width = 14,
  height = 11,
  device = cairo_pdf,
  dpi = 300
)