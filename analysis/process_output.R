library(dplyr)
library(ggpubr)
library(tidyr)
library(forcats)
rm(list=ls())


# Functions ---------------------------------------------------------------

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
             upper = ifelse(metric == "TPR", value, value + diff),
             lower = ifelse(metric == "TPR", value + diff, value),
             norm = fct_rev(norm)
             
      ) |>
      ggplot(aes(y = norm, x = value, color = metric)) +
      geom_vline(xintercept = 0.1, color = "darkred", linetype = "dashed") +
      geom_linerange(
        aes(y = norm, xmax = upper, xmin = lower, group = metric),
        inherit.aes = F,
        color = "gray70",
        linewidth = 0.5
      ) +
      geom_point(aes(x = value), shape = 15, size = 1.2) +
      scale_x_continuous(
        limits = c(0,1),
        breaks = seq(0,1,0.2)
      ) +
      scale_color_manual(
        values = c("red", "blue")
        
      ) +
      labs(x = "Rate",
           color = "Normalization",
           fill = "Normalization",
           title = paste0("n = ", n_use)
           
      ) +
      theme_minimal() +
      facet_grid(rows = vars(variance), cols = vars(prop_signal)) + 
      theme(
        # Panel
        panel.border = element_rect(color = "gray50", fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray90"),
        # panel.grid.minor.y = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        # Text
        plot.title = element_text(size = 11),
        strip.text = element_text(size = 9
        ),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.title = element_blank(),
        axis.title = element_blank(),
        legend.spacing = unit(0.3, "mm"),
        legend.margin = margin(0,0,0,0),
        # Misc
        strip.background = element_rect(fill = "gray90",
                                        color = "gray50"),
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
      mutate(metric = ifelse(metric == "tpr_m", "TPR", "FDR"),
             diff = rep(value[metric == "TPR"] - value[metric != "TPR"], each = 2) * c(-1,1),
             upper = ifelse(metric == "TPR", value, value + diff),
             lower = ifelse(metric == "TPR", value + diff, value),
             norm = fct_rev(norm)
             
      ) |>
      ggplot(aes(y = norm, x = value, color = metric)) +
      geom_vline(xintercept = 0.1, color = "darkred", linetype = "dashed") +
      geom_linerange(
        aes(y = norm, xmax = upper, xmin = lower, group = metric),
        inherit.aes = F,
        color = "gray70",
        linewidth = 0.5
      ) +
      geom_point(aes(x = value), shape = 15, size = 1.2) +
      scale_x_continuous(
        limits = c(0,1),
        breaks = seq(0,1,0.2)
      ) +
      scale_color_manual(
        values = c("red", "blue")
        
      ) +
      labs(x = "Rate",
           color = "Normalization",
           fill = "Normalization",
           title = data_use
      ) +
      theme_minimal() +
      facet_grid(rows = vars(method), cols = vars(prop_signal)) + 
      theme(
        # Panel
        panel.border = element_rect(color = "gray50", fill = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray90"),
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
        legend.box.spacing = unit(0, "pt"),
        legend.margin=margin(c(0,0,0,0)),
        legend.spacing = unit(0, "pt"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        
        # Misc
        strip.background = element_rect(fill = "gray90",
                                        color = "gray50"),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
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




