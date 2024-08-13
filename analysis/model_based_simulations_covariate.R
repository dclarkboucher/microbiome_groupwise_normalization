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
library(MicrobiomeStat)
source("utils/functions.R")


# Simulation --------------------------------------------------------------

# Settings
seed_param <- 1234 # random seed used in parameter generation
ns <- c(200, 500) # sample size
qs <- c(300) # number of taxa
prop_signal <- c(0.1, 0.20, 0.30) # proportion of non-zero log fold changes
alphas <- -qnorm(c(0.5)) # 0.50 = 50% zeros
s2_vs <- c(1, 1.5, 2) # variance parameter
beta1_means <- c(1) # mean of the non-zero log fold changes
methods <- c("edgeR","DESeq2") # DAA methods
norms <- c("FTSS", "G-RLE", "TSS", "GMPR", "Wrench", "TMM", "RLE", "CSS") # norm methods

# Library size parameters
mean_lib <- 60000
var_lib <- 15000^2
size_lib <- mean_lib^2 / (var_lib - mean_lib) # used in NB distribution

# Run simulations
fdr_results <- list()
tpr_results <- list()
j <- 1
for (n in ns){
  
  for(s2_v in s2_vs){
    
    for (q in qs) {
      
      for (ps in prop_signal) {
        
        for (b1_m in beta1_means){
          
          # Sample signals
          set.seed(seed_param)
          beta0 <- rnorm(q, sd = 1)
          q1 <- ceiling(q * ps)
          beta1 <- rep(0, q)
          true_beta <- logical(q)
          if (q1 > 0){
            temp <- rnorm(q1)
            beta1[seq_len(q1)] <- temp + b1_m - mean(temp)
            true_beta[seq_len(q1)] <- TRUE
          }
          
          # Sample covariance matrix
          Sigma_z <-
            cov2cor(LaplacesDemon::rinvwishart(2 * q, diag(q)))
          Sigma_v <- Sigma_z * s2_v
          
          # covariate effects
          beta_c <- ifelse(beta1 == 0, rnorm(q, 0, 0.5), beta1/2)
          
          # Generate X
          x <- rep(1, n)
          x[seq_len(n/2)] <- 0
          
          
          for (a0 in alphas){
            
            for (r in seq_len(total_reps)) {
              
              set.seed(r + 1000)
              
              # Get covariate 
              cov <- x + rnorm(n)
              nu <- exp(beta0 + beta1 %*% t(x) + beta_c %*% t(cov))
              
              # Generate library size
              S <- (rnbinom(n = n, size = size_lib, mu = mean_lib))
              
              # Generate latent variables
              Z <- t(mvtnorm::rmvnorm(n = n, sigma = Sigma_z))
              V <- t(mvtnorm::rmvnorm(n = n, sigma = Sigma_v))
              not_zero <- Z + a0 > 0
              mu <- nu * exp(V) * not_zero
              
              # Sample data
              Y <- matrix(NA, q, n)
              for (i in 1:n) {
                Y[, i] <- as.numeric(rmultinom(1, S[i], prob = mu[, i]))
                
              }
              
              for (norm in norms){
                
                # Calculate normalization factor
                quiet(offset <- get_offset(Y = Y, x = x, method = norm))
                
                for (method in methods){
                  
                  out <- analysis_wrapper_covariate(Y, x, cov, offset, 
                                                    method = method)
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
                      r = r,
                      a0 = a0, 
                      s2_v = s2_v,
                      q = q,
                      q1 = q1,
                      b1_m = b1_m,
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
                  
                  tpr_results[[j]] <-
                    data.frame(
                      method = method,
                      norm = norm,
                      n = n,
                      r = r,
                      a0 = a0, 
                      s2_v = s2_v,
                      q = q,
                      q1 = q1,
                      b1_m = b1_m,
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
    }
  }
}

fdr_results_sum <- bind_rows(fdr_results)
tpr_results_sum <- bind_rows(tpr_results)
results_sum <- list(fdr = fdr_results_sum,
                    tpr = tpr_results_sum)

filename <- paste0("output/model_covariate_sims.rda")
save(results_sum, file = filename)
