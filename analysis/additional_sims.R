library(tidyr)
library(tibble)
library(dplyr)
library(edgeR)
library(DESeq2)
library(metagenomeSeq)
library(MicrobiomeStat)
source("utils/general_sim_functions.R")

total_reps <- 4
reps_do <- seq_len(total_reps)


# Null signal sims --------------------------------------------------------
null_sims_data <-
  expand_grid(
    data = c("HPFS","PHACS"),
    ps = 0,
    # norm = c("TSS","RLE", "GMPR", "Wrench",
    #          "G-RLE", "FTSS"),
    norm = c("Wrench",
             "G-RLE", "FTSS"),
    method = c("mgs", "DESeq2"),
    signal = "right" # N(1,1) signals
  )

out <-
  general_data_sims(
    settings = null_sims_data,
    reps_do = reps_do
)

save(out, file = paste0("output/null_signal_sims.rda"))


# Little comp. bias sims --------------------------------------------------
less_bias_sims_data <-
  expand_grid(
    data = c("HPFS","PHACS"),
    ps = c(0.2),
    norm = c("TSS","RLE", "GMPR", "Wrench",
             "G-RLE", "FTSS"),
    method = c("mgs", "DESeq2"),
    signal = c("mixed", "center", "right")
  )

out <-
  general_data_sims(
    settings = less_bias_sims_data,
    reps_do = reps_do
  )

save(out, file = paste0("output/less_bias_sims.rda"))




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

save(out, file = paste0("output/libsize_sims.rda"))


# LinDA sims --------------------------------------------------------------

# Comparison to DAA method LinDA
linda_sims <-
  expand_grid(
    n = c(200, 500),
    ps = c(0.1, 0.2, 0.3),
    s2_v = 2,
    norm = c("linda", "GMPR", "Wrench", "G-RLE", "FTSS"),
    # Note: if LinDA is provided as a norm method to this function, 
    # the DAA method (e.g., mgs) is ignored. 
    method = c("mgs"),
    signal = c("right")
  )


out <-
  general_model_sims(
    settings = linda_sims,
    reps_do = reps_do
  )

save(out, file = paste0("output/linda_sims.rda"))

