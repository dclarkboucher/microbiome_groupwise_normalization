# script

# other packages: LaplacesDemon
rm(list=ls())
total_reps <- 1 

# Replicate the analysis by setting total_reps = 1000 and splitting the 
# simulations into separate, parallel jobs

# Packages ----------------------------------------------------------------

library(tibble)
library(dplyr)
library(edgeR)
library(DESeq2)
library(metagenomeSeq)

# Generate toy data -------------------------------------------------------
n <- 200
q <- 300
beta0 <- rnorm(q)
s2_v <- 1.5
Sigma_z <- cov2cor(LaplacesDemon::rinvwishart(2 * q, diag(q)))
Sigma_v <- Sigma_z * s2_v
Z <- t(mvtnorm::rmvnorm(n = n, sigma = Sigma_z))
V <- t(mvtnorm::rmvnorm(n = n, sigma = Sigma_v))
mu <- exp(beta0 + V) * (Z > 0)
mean_lib <- 60000
var_lib <- 15000^2
size_lib <- mean_lib^2 / (var_lib - mean_lib)
S <- rnbinom(n = n, size = size_lib, mu = mean_lib)
Y <- matrix(NA, nrow = q, ncol = n)
for (i in seq_len(n)) Y[, i] <- as.numeric(rmultinom(1, S[i], prob = mu[, i]))
dataset_lists <- list()
dataset_lists[["toy_data"]] <- list(Y = Y, S = S)
rm(list = setdiff(ls(), c("dataset_lists", "total_reps")))


# Simulation --------------------------------------------------------------

# Load some functions
source("utils/functions.R")

# Settings
seed_param <- 1234 # random seed used in parameter generation
prop_signal <- c(0.1, 0.20, 0.30) # proportion of non-zero log fold changes
beta1_means <- c(1) # mean of the non-zero log fold changes
methods <- c("edgeR","metagenomeSeq","DESeq2") # DAA methods
norms <- c("FTSS", "G-RLE", "TSS", "GMPR", "Wrench", "TMM", "RLE", "CSS") # norm methods
datasets <- c("toy_data") # c("PHACS", "MLVS/MBS")

# Run simulations
fdr_results <- list()
fpr_results <- list()
j <- 1

for (d in datasets){
  
  set.seed(1)
  Y_true <- dataset_lists[[d]]$Y
  S <- dataset_lists[[d]]$S
  q <- nrow(Y_true)
  n <- ncol(Y_true)
  Y_true <- Y_true[sample(q),]
  
  for (ps in prop_signal) {
    
    for (b1_m in beta1_means){
      
      # Sample signals
      set.seed(seed_param)
      beta0 <- rnorm(q, sd = 1) # don't use these
      q1 <- ceiling(q * ps)
      beta1 <- rep(0, q)
      true_beta <- logical(q)
      if (q1 > 0){
        temp <- rnorm(q1)
        beta1[seq_len(q1)] <- temp + b1_m - mean(temp)
        true_beta[seq_len(q1)] <- TRUE
      }
      
     for (r in seq_len(total_reps)) {
        
        set.seed(r + 1000)
        
        # Generate Data
        synth_data <- get_synthetic_data(Y = Y_true, beta1 = beta1, S = S)
        Y <- synth_data$Y
        x <- synth_data$x
        
        for (norm in norms){
          
          # Calculate normalization factor
          quiet(offset <- get_offset(Y = Y, x = x, method = norm))
          
          for (method in methods){
            
            out <- analysis_wrapper(Y, x, offset, method = method)
            beta1_hat <- out$beta1_hat
            pv <- out$pv
            
            # FDR results 
            cutoffs <- seq(0.01, 0.2, 0.01) 
            nc <- length(cutoffs)
            fdp <- fpr <- tpr <- numeric(nc)
            pv_adj <- p.adjust(pv, method = "BH")
            for (l in seq_len(nc)){
              cut <- cutoffs[l]
              positive <- pv_adj < cut
              true_positive <- positive & true_beta
              false_positive <- positive & !true_beta
              
              tpr[l] <- sum(true_positive) / q1
              fpr[l] <- sum(false_positive) / (q - q1)
              
              if (!any(positive)) {
                fdp[l] <- 0
              } else {
                fdp[l] <- sum(positive & (!true_beta)) / sum(positive)
                
              }
            }
            
            fdr_results[[j]] <-
              data.frame(
                method = method,
                norm = norm,
                n = n,
                data = d,
                q = q,
                q1 = q1,
                b1_m = b1_m,
                r = r,
                fdr = cutoffs, # true FDR
                tpr = tpr, 
                fpr = fpr,
                fdp = fdp # observed false discovery proportion (FDP)
              )
            
            # TPR results
            
            false_pv_indices <- which(!(true_beta[order(pv)]))
            fp <- seq_len(q - q1)
            fpr <- fp / (q- q1)
            tp <- false_pv_indices - fp
            tpr <- tp / q1
            
            fpr_results[[j]] <-
              data.frame(
                method = method,
                norm = norm,
                n = n,
                data = d,
                q = q,
                q1 = q1,
                b1_m = b1_m,
                r = r,
                tpr = tpr, 
                fpr = fpr,
                tp = tp,
                fp = fp
              )
            
            j <- j + 1
          }
        }
      }
      
    }
  }
  
}

fdr_results_sum <- bind_rows(fdr_results)
fpr_results_sum <- bind_rows(fpr_results)
results_sum <- list(fdr = fdr_results_sum,
                    fpr = fpr_results_sum)

filename <- paste0("output/synthetic_data_sims.rda")
save(results_sum, file = filename)
