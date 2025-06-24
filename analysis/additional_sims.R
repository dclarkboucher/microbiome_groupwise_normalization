library(tidyr)
library(tibble)
library(dplyr)
library(edgeR)
library(DESeq2)
library(metagenomeSeq)
library(MicrobiomeStat)
source("utils/general_sim_functions.R")

total_reps <- 1
njobs <- 50

# Null signal sims --------------------------------------------------------
null_sims_data <-
  expand_grid(
    data = c("HPFS","PHACS"),
    ps = 0,
    norm = c("TSS","RLE", "GMPR", "Wrench",
             "G-RLE", "FTSS"),
    method = c("mgs", "DESeq2"),
    signal = "right"
  )

out1 <-
  general_data_sims(
    settings = null_sims_data,
    reps_do = reps_do
)

save(out, file = paste0("~/output/null_signal_sims.rda"))


# Little comp. bias sims --------------------------------------------------
signal_sims_data <-
  expand_grid(
    data = c("HPFS","PHACS"),
    ps = c(0.2, 0.3),
    norm = c("TSS","RLE", "GMPR", "Wrench",
             "G-RLE", "FTSS"),
    method = c("mgs", "DESeq2"),
    signal = c("mixed", "center", "right")
  )

out <-
  general_data_sims(
    settings = signal_sims_data,
    reps_do = reps_do
  )

save(out, file = paste0("~/output/less_bias_sims.rda"))




# Imbalanced library size -------------------------------------------------

# Library sizes in the experimental group (x=1) are multiplied by a scalar
libsize_sims_model <-
  expand_grid(
    n = c(200, 500),
    ps = c(0.2),
    s2_v = 2,
    norm = c("TSS","RLE", "GMPR", "Wrench",
             "G-RLE", "FTSS"),
    method = c("mgs", "DESeq2"),
    signal = c("right"),
    libscale = c(.90, 1, 1.10) # library size scalar
  )


out <-
  general_model_sims(
    settings = libsize_sims_model,
    reps_do = reps_do
  )

save(out, file = paste0("~/output/libsize_sims.rda"))




