library(dplyr)
library(ggpubr)
library(tidyr)
library(forcats)
rm(list=ls())


# Functions ---------------------------------------------------------------
largefont <- 7
smallfont <- 6
plot_method_n <- function(method_use, n_use, type = c("fdr", "fpr"), results){
  
  type <- match.arg(type)
  
  if (type == "fpr"){
    
    results$fpr |>
      filter(
        method == method_use,
        n == n_use
      ) |>
      ggplot(aes(x = fpr, y = tpr_m,
                 color = factor(norm, ordered = F)
      )) +
      geom_line(linewidth = 0.5) +
      scale_y_continuous(
        limits = c(0,1),
        breaks = seq(0,1,0.2)
      ) +
      scale_x_continuous(
        limits = c(0,0.4),
        breaks = seq(0,0.4,0.1)
      ) +
      labs(x = "Fixed FPR", y = "Mean TPR",
           color = "Normalization",
           fill = "Normalization",
           title = paste0("n = ", n_use)
      ) +
      theme_minimal() +
      facet_grid(rows = vars(variance), 
                 cols = vars(prop_signal),
                 margins = F) + 
      theme(
        # Panel
        panel.border = element_rect(color = "gray50", fill = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80",
                                          linewidth = 0.3),
        # panel.grid.minor.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        # Text
        plot.title = element_text(size = 11),
        strip.text = element_text(size = 9
        ),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.margin = margin(0,0,0,0), 
        legend.box.spacing = unit(0.5, "mm"),
        
        # Misc
        strip.background = element_rect(fill = NA,
                                        color = NA),
        axis.ticks = element_blank(),
        legend.position = "right"
      )
    
    
  } else if (type == "fdr"){
    
    results$fdr |>
      filter(n == n_use,
             method == method_use,
             fdr < 0.11 & fdr > 0.09
      ) |>
      pivot_longer(cols = c("tpr_m", "fdp_m"),
                   names_to = "metric",
                   values_to = "value"
      ) |>
      mutate(metric = ifelse(metric == "tpr_m", "TPR", "FDR"),
             diff = rep(value[metric == "TPR"] - value[metric != "TPR"], each = 2) * c(-1,1),
             upper = value + 2 * ifelse(metric == "TPR", tpr_se, fdp_se),
             lower = value - 2 * ifelse(metric == "TPR", tpr_se, fdp_se),
             norm_lab = fct_rev(norm_lab)
             
      ) |>
      ggplot(aes(y = norm_lab, x = value, color = metric)) +
      geom_vline(xintercept = fdr_cutoff, color = "black") +
      geom_col(
        aes(y = norm_lab, x = value, group = metric, fill = metric),
        position = position_dodge(0.6),
        width = 0.6, color = NA
      ) +
      geom_errorbar(
        aes(y = norm_lab, xmax = upper, xmin = lower, group = metric),
        position = position_dodge(width = 0.6),
        color = "gray40",
        width = 0.3
      ) +
      scale_x_continuous(
        limits = c(0,1),
        breaks = seq(0,0.8,0.2),
        expand = c(0,0),
        labels = paste0(round(seq(0,0.8,0.2),1))
      ) +
      scale_y_discrete(expand = expansion(add = c(0.4,0.4))) + 
      scale_fill_manual(
        values = c("#e52b50", "#2bade5")
        
      ) +
      labs(x = "Rate",
           color = "Normalization",
           fill = "Normalization"
           
      ) +
      theme_minimal() +
      facet_grid(cols = vars(variance1), 
                 rows = vars(prop_signal),
                 margins = F) + 
      theme(
        # Panel
        panel.border = element_rect(color = "gray50", fill = NA),
        axis.ticks.x = element_line(color = "gray50", linewidth = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "solid"),
        # panel.grid.minor.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        # Text
        # plot.title = element_text(size = title_font, hjust = 0.5),
        strip.text = element_text(size = largefont
        ),
        axis.text = element_text(size = smallfont),
        legend.text = element_text(size = smallfont),
        legend.title = element_blank(),
        axis.title = element_blank(),
        legend.spacing = unit(0, "mm"),
        legend.key.width = unit(0,"mm"),
        legend.margin = margin(0,0,0,0),
        plot.margin = margin(1,1,1,1),
        # Misc
        strip.background = element_rect(fill = NA,
                                        color = NA),
        axis.ticks.y = element_blank(),
        legend.position = "right"
      )
    
    
  }
  
}


plot_method <- function(method, type = c("fdr", "fpr"), results){
  
  type <- match.arg(type)
  p1 <- plot_method_n(method, 200, type = type, results = results)
  p2 <- plot_method_n(method, 500, type = type, results = results)
  ggarrange(p1, p2, ncol = 1, common.legend = TRUE,
            legend = "bottom")
  
}

plot_synthetic_data <- function(data_use, type = c("fdr", "fpr"), results){
  
  if (type == "fpr"){
    results$fpr |>
      filter(data == data_use) |>
      ggplot(aes(x = fpr, y = tpr_m,
                 color = factor(norm, ordered = F)
      )) +
      geom_line(linewidth = 0.5) +
      scale_y_continuous(
        limits = c(0,1),
        breaks = seq(0,1,0.2)
      ) +
      scale_x_continuous(
        limits = c(0,0.4),
        breaks = seq(0,0.4,0.1)
      ) +
      labs(x = "Fixed FPR", y = "Mean TPR",
           color = "Normalization",
           fill = "Normalization",
           title = data_use
      ) +
      theme_minimal() +
      facet_grid(rows = vars(method), 
                 cols = vars(prop_signal),
                 margins = F) + 
      theme(
        # Panel
        panel.border = element_rect(color = "gray50", fill = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80",
                                          linewidth = 0.3),
        # panel.grid.minor.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        # Text
        plot.title = element_text(size = 11),
        strip.text = element_text(size = 9
        ),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.margin = margin(0,0,0,0), 
        legend.box.spacing = unit(0.5, "mm"),
        
        # Misc
        strip.background = element_rect(fill = NA,
                                        color = NA),
        axis.ticks = element_blank(),
        legend.position = "right"
      )
  } else if (type == "fdr"){
    
    fdr |>
      filter(data == data_use,
             fdr < 0.11 & fdr > 0.09
      ) |>
      pivot_longer(cols = c("tpr_m", "fdp_m"),
                   names_to = "metric",
                   values_to = "value"
      ) |>
      mutate(
        metric = ifelse(metric == "tpr_m", "TPR", "FDR"),
        
        upper = ifelse(metric == "TPR", value + 2 * fdp_se, value + 2 * tpr_se),
        lower = ifelse(metric == "TPR", value - 2 * fdp_se, value + 2 * fdp_se),
        norm_lab = fct_rev(norm_lab)
        
      ) |>
      ggplot(aes(y = norm_lab, x = value, color = metric)) +
      geom_vline(xintercept = fdr_cutoff, color = "black") +
      geom_col(
        aes(y = norm_lab, x = value, group = metric, fill = metric),
        position = position_dodge(0.6),
        width = 0.6, color = NA
      ) +
      geom_errorbar(
        aes(y = norm_lab, xmax = upper, xmin = lower, group = metric),
        position = position_dodge(width = 0.6),
        color = "gray40",
        width = 0.3
      ) + 
      # geom_point(aes(x = value), shape = 15, size = 1.2) +
      scale_x_continuous(
        limits = c(0,1),
        breaks = seq(0,0.8,0.2),
        expand = c(0,0),
        labels = paste0(round(seq(0,0.8,0.2),1))
      ) +
      scale_y_discrete(expand = expansion(add = c(0.4,0.4))) + 
      scale_fill_manual(
        values = c("#e52b50", "#2bade5")
        
      ) +
      labs(x = "Rate",
           color = "Normalization",
           fill = "Normalization"
           
      ) +
      theme_minimal() +
      facet_grid(cols = vars(d, method), 
                 rows = vars(prop_signal),
                 margins = F) + 
      theme(
        # Panel
        panel.border = element_rect(color = "gray50", fill = NA),
        axis.ticks.x = element_line(color = "gray50", linewidth = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "solid"),
        # panel.grid.minor.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        # Text
        # plot.title = element_text(size = title_font, hjust = 0.5),
        strip.text = element_text(size = largefont
        ),
        axis.text = element_text(size = smallfont),
        legend.text = element_text(size = smallfont),
        legend.title = element_blank(),
        axis.title = element_blank(),
        legend.spacing = unit(0, "mm"),
        legend.key.width = unit(0,"mm"),
        legend.margin = margin(0,0,0,0),
        plot.margin = margin(1,1,1,1),
        # Misc
        strip.background = element_rect(fill = NA,
                                        color = NA),
        axis.ticks.y = element_blank(),
        legend.position = "right"
      )
    
    
    
  }
}


# Model-based simulations -------------------------------------------------
load("output/model_sims.rda")


# Process FDR
fdr <- 
  results_sum$fdr |> 
  group_by(method, norm, n, a0, s2_v, q, q1, b1_m, fdr) |> 
  summarize(across(c(tpr, fpr, fdp), 
                   .fns = 
                     list(
                       m = mean, 
                       med = median, 
                       se = ~ sd(.x) / sqrt(length(.x)),
                       l = ~ quantile(.x, 0.025),
                       u = ~ quantile(.x, 0.975)                  )
  )) |>
  ungroup() |> 
  mutate(
    prop_signal = paste0(floor(q1 / q * 100), "% signals"),
    variance = factor(as.integer(factor(s2_v)), labels = c("Low Var.", "Medium Var.", "High Var.")),
    norm = factor(norm, 
                  levels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  labels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  ordered = TRUE
    )
    
  )

fpr <- 
  results_sum$fpr |> 
  group_by(method, norm, n, a0, s2_v, q, q1, b1_m, fpr) |> 
  summarize(across(c(tpr), 
                   .fns = 
                     list(
                       m = mean, 
                       med = median, 
                       se = ~ sd(.x) / sqrt(length(.x)),
                       l = ~ quantile(.x, 0.025),
                       u = ~ quantile(.x, 0.975))
  )) |>
  ungroup() |> 
  mutate(
    prop_signal = paste0(floor(q1 / q * 100), "% signals"),
    variance = factor(as.integer(factor(s2_v)), labels = c("Low Var.", "Medium Var.", "High Var.")),
    norm = factor(norm, 
                  levels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  labels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  ordered = TRUE
    )
    
  )

# Plot results
dat <- list(fdr = fdr, fpr = fpr)
for (method in c("edgeR", "DESeq2", "metagenomeSeq")){
  
  p <- plot_method(method, "fdr", dat)
  ggsave(
    filename = paste0("figures/model_sims_fdr_",method,".jpeg"),
    plot = p,
    width = 6,
    height = 9, 
    units = "in"
  )
  
  p <- plot_method(method, "fpr", dat)
  ggsave(
    filename = paste0("figures/model_sims_fpr_",method,".jpeg"),
    plot = p,
    width = 6,
    height = 9, 
    units = "in"
  )
  
}






# Model-based simulations (confounder) ------------------------------------

load("output/model_confounder_sims.rda")

# Process FDR
fdr <- 
  results_sum$fdr |> 
  group_by(method, norm, n, a0, s2_v, q, q1, b1_m, fdr) |> 
  summarize(across(c(tpr, fpr, fdp), 
                   .fns = 
                     list(
                       m = mean, 
                       med = median, 
                       se = ~ sd(.x) / sqrt(length(.x)),
                       l = ~ quantile(.x, 0.025),
                       u = ~ quantile(.x, 0.975)                  )
  )) |>
  ungroup() |> 
  mutate(
    prop_signal = paste0(floor(q1 / q * 100), "% signals"),
    variance = factor(as.integer(factor(s2_v)), labels = c("Low Var.", "Medium Var.", "High Var.")),
    norm = factor(norm, 
                  levels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  labels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  ordered = TRUE
    )
    
  )

fpr <- 
  results_sum$fpr |> 
  group_by(method, norm, n, a0, s2_v, q, q1, b1_m, fpr) |> 
  summarize(across(c(tpr), 
                   .fns = 
                     list(
                       m = mean, 
                       med = median, 
                       se = ~ sd(.x) / sqrt(length(.x)),
                       l = ~ quantile(.x, 0.025),
                       u = ~ quantile(.x, 0.975))
  )) |>
  ungroup() |> 
  mutate(
    prop_signal = paste0(floor(q1 / q * 100), "% signals"),
    variance = factor(as.integer(factor(s2_v)), labels = c("Low Var.", "Medium Var.", "High Var.")),
    norm = factor(norm, 
                  levels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  labels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  ordered = TRUE
    )
    
  )

# Plot results
dat <- list(fdr = fdr, fpr = fpr)
for (method in c("edgeR", "DESeq2", "metagenomeSeq")){
  
  p <- plot_method(method, "fdr", dat)
  ggsave(
    filename = paste0("figures/model_confounder_sims_fdr_",method,".jpeg"),
    plot = p,
    width = 6,
    height = 9, 
    units = "in"
  )
  
  p <- plot_method(method, "fpr", dat)
  ggsave(
    filename = paste0("figures/model_confounder_fpr_",method,".jpeg"),
    plot = p,
    width = 6,
    height = 9, 
    units = "in"
  )
  
}


# Synthetic data simulations ----------------------------------------------

load("output/synthetic_data_sims.rda")

fdr <- 
  results_sum$fdr |> 
  group_by(method, norm, n, data, q, q1, b1_m, fdr) |> 
  summarize(across(c(tpr, fpr, fdp), 
                   .fns = 
                     list(
                       m = mean, 
                       med = median, 
                       se = ~ sd(.x) / sqrt(length(.x)),
                       l = ~ quantile(.x, 0.025),
                       u = ~ quantile(.x, 0.975)                  )
  ))  |>
  ungroup() |> 
  mutate(
    prop_signal = paste0(floor(q1 / q * 100), "% signals"),
    norm = factor(norm, 
                  levels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  labels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  ordered = TRUE
    )
    
  )

fpr <- 
  results_sum$fpr |> 
  group_by(method, norm, n, data, q, q1, b1_m, fpr) |> 
  summarize(across(tpr, 
                   .fns = 
                     list(
                       m = mean, 
                       med = median, 
                       se = ~ sd(.x) / sqrt(length(.x)),
                       l = ~ quantile(.x, 0.025),
                       u = ~ quantile(.x, 0.975)                  )
  ))  |>
  ungroup() |> 
  mutate(
    prop_signal = paste0(floor(q1 / q * 100), "% signals"),
    norm = factor(norm, 
                  levels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  labels = c("FTSS","G-RLE","GMPR", "Wrench", "RLE","CSS", "TMM", "TSS"),
                  ordered = TRUE
    )
    
  )

# Plot results for toy dataset
results <- list(fpr = fpr, fdr = fdr)

data_use <- "toy_data"
p <- plot_synthetic_data(data_use = data_use, type = "fpr", results = results)
ggsave(
  filename = paste0("figures/synthetic_fpr_",data_use,".jpeg"),
  plot = p,
  width = 6,
  height = 9, 
  units = "in"
)

p <- plot_synthetic_data(data_use = data_use, type = "fdr", results = results)
ggsave(
  filename = paste0("figures/synthetic_fdr_",data_use,".jpeg"),
  plot = p,
  width = 6,
  height = 9, 
  units = "in"
)






# Global null simulations -------------------------------------------------

load("output/null_signal_sims.rda")
out$fdr |> 
  filter(fdr == 0.05) |> 
  mutate(fp = q * fdp) |> # fdp: observed false positive proportion 
  group_by(method, norm, data, q) |> 
  summarize(
    fpr_m = mean(fpr),
    fpr_se = sd(fpr) / sqrt(n())
  ) |> 
  mutate(
    norm = fct_rev(factor(norm, 
                          levels = c("FTSS", "G-RLE", "GMPR", 
                                     "Wrench", "RLE", "TSS"),
                          labels = c("FTSS", "G-RLE", "GMPR", 
                                     "Wrench", "RLE", "TSS"))),
    data = ifelse(data == "HPFS", "MLVS/MBS", data),
    method = ifelse(method == "mgs", "MetagenomeSeq", method),
    fp = fpr_m * q,
    fp_lower = (fpr_m - 2 * fpr_se) * q, 
    fp_upper  = (fpr_m + 2 * fpr_se) * q
  ) |> 
  select(method, norm, data, fp, fp_lower, fp_upper) |> 
  ggplot(aes(y = norm, x = fp, fill = method)) +
  facet_wrap(~data, scales = "free_x") +
  geom_col(position = position_dodge(width = 0.6),
           width = 0.6
  ) + 
  geom_errorbar(
    aes(y = norm, xmax = fp_upper,
        xmin = fp_lower, group = method),
    position = position_dodge(width = 0.6),
    color = "gray30",
    width = 0.4,
    linewidth = 0.4
  ) + 
  scale_fill_manual(
    values = c("#e52b50", "#2bade5")
    
  ) +
  theme_minimal() + 
  labs(x = "False positives") + 
  theme(
    # Panel
    panel.border = element_rect(color = "gray50", fill = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    
    # Text
    strip.text = element_text(size = largefont),
    axis.text = element_text(size = smallfont),
    legend.text = element_text(size = smallfont),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = largefont),
    legend.spacing = unit(0, "mm"),
    legend.key.width = unit(5,"mm"),
    legend.key.height = unit(4,"mm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(1,1,1,1),
    # Misc
    strip.background = element_rect(fill = NA,
                                    color = NA),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )


# Imbalanced library size simulations -------------------------------------

load("output/libsize_sims.rda")
out$fdr |> 
  group_by(method, norm, n, signal, libscale, s2_v, q1, q, fdr) %>% 
  summarize(
    fdp_m = mean(fdp),
    fdp_se = sd(fdp) / sqrt(n()),
    fdp_l = as.numeric(quantile(fdp, 0.025)),
    fdp_u = as.numeric(quantile(fdp, 0.975)),
    fpr_m = mean(fpr),
    fpr_se = sd(fpr) / sqrt(n()),
    fpr_l = as.numeric(quantile(fpr, 0.025)),
    fpr_u = as.numeric(quantile(fpr, 0.975)),
    tpr_m = mean(tpr),
    tpr_se = sd(tpr) / sqrt(n()),
    tpr_l = as.numeric(quantile(tpr, 0.025)),
    tpr_u = as.numeric(quantile(tpr, 0.975))
  ) |> 
  filter(fdr == 0.05, libscale != 1) |> 
  mutate(
    norm = fct_rev(factor(norm, 
                          levels = c("FTSS", "G-RLE", "GMPR", 
                                     "Wrench", "RLE", "TSS"),
                          labels = c("FTSS", "G-RLE", "GMPR", 
                                     "Wrench", "RLE", "TSS"))),
    tpr_lower = tpr_m - 2 * tpr_se,
    tpr_upper = tpr_m + 2 * tpr_se,
    fdp_lower = fdp_m - 2 * fdp_se,
    fdp_upper = fdp_m + 2 * fdp_se,
    method = ifelse(method == "mgs", "MetagenomeSeq", "DESeq2"),
    libscale = paste0("Library scalar: ", libscale),
    n = paste0("n = ",n)
  ) |> 
  select(method, norm, n, libscale, fdp_m, fdp_lower, fdp_upper, 
         tpr_m, tpr_lower, tpr_upper) |> 
  pivot_longer(
    cols = c("fdp_m", "tpr_m"),
    names_to = "metric", 
    values_to = "performance"
  ) |> 
  mutate(
    lower = ifelse(metric == "fdp_m", fdp_lower, tpr_lower),
    upper = ifelse(metric == "fdp_m", fdp_upper, tpr_upper),
    metric = ifelse(metric == "fdp_m", "FDR", "TPR")
  ) |> 
  select(-starts_with("fdp"), -starts_with("tpr")) |> 
  ggplot(aes(y = norm, x = performance, fill = metric)) +
  geom_vline(xintercept = 0.05, color = "black") +
  geom_col(position = position_dodge(0.5),
           width = 0.5) + 
  scale_fill_manual(
    values = c("#e52b50", "#2bade5")
    
  ) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(0.5),
                width = 0.3
  ) + 
  facet_grid(
    cols = vars(n, method),
    rows = vars(libscale)
  ) + 
  theme_minimal() + 
  theme(
    # Panel
    panel.border = element_rect(color = "gray50", fill = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    
    # Text
    strip.text = element_text(size = largefont),
    axis.text = element_text(size = smallfont),
    legend.text = element_text(size = smallfont),
    legend.title = element_blank(),
    axis.title = element_blank(),
    legend.spacing = unit(0, "mm"),
    # legend.key.width = unit(0,"mm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(1,1,1,1),
    # Misc
    strip.background = element_rect(fill = NA,
                                    color = NA),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

# Little comp. bias sims --------------------------------------------------

load("output/less_bias_sims.rda")
out$fdr |> 
  group_by(method, norm, signal, data, q1, q, fdr) %>% 
  summarize(
    fdp_m = mean(fdp),
    fdp_se = sd(fdp) / sqrt(n()),
    fdp_l = as.numeric(quantile(fdp, 0.025)),
    fdp_u = as.numeric(quantile(fdp, 0.975)),
    fpr_m = mean(fpr),
    fpr_se = sd(fpr) / sqrt(n()),
    fpr_l = as.numeric(quantile(fpr, 0.025)),
    fpr_u = as.numeric(quantile(fpr, 0.975)),
    tpr_m = mean(tpr),
    tpr_se = sd(tpr) / sqrt(n()),
    tpr_l = as.numeric(quantile(tpr, 0.025)),
    tpr_u = as.numeric(quantile(tpr, 0.975))
  ) |> 
  filter(
    fdr == 0.05, 
    signal != "right", 
    round(q1/q,2) == 0.2
  ) |> 
  mutate(
    norm = fct_rev(factor(norm, 
                          levels = c("FTSS", "G-RLE", "GMPR", 
                                     "Wrench", "RLE", "TSS"),
                          labels = c("FTSS", "G-RLE", "GMPR", 
                                     "Wrench", "RLE", "TSS"))),
    tpr_lower = tpr_m - 2 * tpr_se,
    tpr_upper = tpr_m + 2 * tpr_se,
    fdp_lower = fdp_m - 2 * fdp_se,
    fdp_upper = fdp_m + 2 * fdp_se,
    prop_signal = paste0(floor(q1/q * 100),"% signals"),
    signal = ifelse(signal == "center", "N(0,1) Signals", "N(-1,1) or N(1,1) Signals"),
    method = ifelse(method == "mgs", "MetagenomeSeq", "DESeq2"),
    method_data = paste0(method," on ",data)
  ) |> 
  select(
    method_data, norm, signal, prop_signal, fdp_m, fdp_lower, fdp_upper, 
    tpr_m, tpr_lower, tpr_upper
  ) |> 
  pivot_longer(
    cols = c("fdp_m", "tpr_m"),
    names_to = "metric", 
    values_to = "performance"
  ) |> 
  mutate(
    lower = ifelse(metric == "fdp_m", fdp_lower, tpr_lower),
    upper = ifelse(metric == "fdp_m", fdp_upper, tpr_upper),
    metric = ifelse(grepl("tpr", metric), "TPR", "FDR")
  ) |> 
  select(-starts_with("fdp"), -starts_with("tpr")) |> 
  ggplot(aes(y = norm, x = performance, fill = metric)) +
  geom_vline(xintercept = 0.05, color = "black") +
  geom_col(position = position_dodge(0.6),
           width = 0.6) + 
  geom_errorbar(
    aes(y = norm, xmax = upper, 
        xmin = lower),
    position = position_dodge(width = 0.6),
    color = "gray30",
    width = 0.4,
    linewidth = 0.4
  ) + 
  scale_fill_manual(
    values = c("#e52b50", "#2bade5")
    
  ) +
  facet_grid(
    cols = vars(method_data),
    rows = vars(signal)
  ) + 
  theme_minimal() + 
  theme(
    # Panel
    panel.border = element_rect(color = "gray50", fill = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    
    # Text
    strip.text = element_text(size = largefont),
    axis.text = element_text(size = smallfont),
    legend.text = element_text(size = smallfont),
    legend.title = element_blank(),
    axis.title = element_blank(),
    legend.spacing = unit(0, "mm"),
    # legend.key.width = unit(0,"mm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(1,1,1,1),
    # Misc
    strip.background = element_rect(fill = NA,
                                    color = NA),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )



# LinDA sims --------------------------------------------------------------

load("output/linda_sims.rda")
out$fdr |> 
  group_by(method, norm, n, signal, libscale, s2_v, q1, q, fdr) %>% 
  summarize(
    fdp_m = mean(fdp),
    fdp_se = sd(fdp) / sqrt(n()),
    fdp_l = as.numeric(quantile(fdp, 0.025)),
    fdp_u = as.numeric(quantile(fdp, 0.975)),
    fpr_m = mean(fpr),
    fpr_se = sd(fpr) / sqrt(n()),
    fpr_l = as.numeric(quantile(fpr, 0.025)),
    fpr_u = as.numeric(quantile(fpr, 0.975)),
    tpr_m = mean(tpr),
    tpr_se = sd(tpr) / sqrt(n()),
    tpr_l = as.numeric(quantile(tpr, 0.025)),
    tpr_u = as.numeric(quantile(tpr, 0.975))
  ) |> 
  ungroup() |> 
  filter(fdr == 0.05) |> 
  mutate(
    prop_signal = paste0(round(100 * q1/q),"% signals"),
    norm_lab = ifelse(norm == "linda", "LinDA", paste("MetagenomeSeq:", norm)),
    tpr_upper = tpr_m + 2 * tpr_se,
    tpr_lower = tpr_m - 2 * tpr_se,
    fdp_lower = fdp_m - 2 * fdp_se,
    fdp_upper = fdp_m + 2 * fdp_se,
    method = ifelse(method == "mgs", "MetagenomeSeq", "DESeq2"),
    n = paste0("n = ",n)
  ) |> 
  select(norm_lab, n, prop_signal, fdp_m, fdp_lower, fdp_upper, 
         tpr_m, tpr_lower, tpr_upper) |> 
  pivot_longer(
    cols = c("fdp_m", "tpr_m"),
    names_to = "metric"
  ) |> 
  mutate(
    lower = ifelse(metric == "fdp_m", fdp_lower, tpr_lower),
    upper = ifelse(metric == "fdp_m", fdp_upper, tpr_upper),
    metric = ifelse(metric == "fdp_m", "FDR", "TPR")
  ) |> 
  ggplot(aes(y = norm_lab, x = value, color = metric)) +
  geom_vline(xintercept = 0.05, color = "black") +
  geom_col(
    aes(y = norm_lab, x = value, group = metric, fill = metric),
    position = position_dodge(0.6),
    width = 0.6, color = NA
  ) +
  geom_errorbar(
    aes(y = norm_lab, xmax = upper, xmin = lower, group = metric),
    position = position_dodge(width = 0.6),
    color = "gray40",
    width = 0.3
  ) +
  scale_x_continuous(
    limits = c(0,1),
    breaks = seq(0,0.8,0.2),
    expand = c(0,0),
    labels = paste0(round(seq(0,0.8,0.2),1))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.4,0.4))) + 
  scale_fill_manual(
    values = c("#e52b50", "#2bade5")
  ) +
  labs(x = "Rate",
       color = "Normalization",
       fill = "Normalization"
       
  ) +
  theme_minimal() +
  facet_grid(rows = vars(n), 
             cols = vars(prop_signal),
             margins = F) + 
  theme(
    # Panel
    panel.border = element_rect(color = "gray50", fill = NA),
    axis.ticks.x = element_line(color = "gray50", linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linetype = "solid"),
    # panel.grid.minor.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # Text
    # plot.title = element_text(size = title_font, hjust = 0.5),
    strip.text = element_text(size = largefont
    ),
    axis.text = element_text(size = smallfont),
    legend.text = element_text(size = smallfont),
    legend.title = element_blank(),
    axis.title = element_blank(),
    legend.spacing = unit(0, "mm"),
    legend.key.width = unit(0,"mm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(1,1,1,1),
    # Misc
    strip.background = element_rect(fill = NA,
                                    color = NA),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )




