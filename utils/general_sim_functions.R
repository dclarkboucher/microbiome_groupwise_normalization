library(tibble)
library(dplyr)
library(edgeR)
library(DESeq2)
library(metagenomeSeq)
library(MicrobiomeStat)

# Helper functions --------------------------------------------------------

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

gmpr <- function (comm, intersect.no = 10, ct.min = 1, trace = TRUE) {
  # Computes the GMPR size factor
  #
  # Args:
  #   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
  #   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
  #   ct.min: the minimum number of counts required to calculate ratios
  
  #
  # Returns:
  #   a vector of the size factors with attribute 'NSS'. Samples with distinct sets of features will be output as NA.
  #         NSS:   number of samples with significant sharing (> intersect.no) including itself
  
  # mask counts < ct.min
  comm[comm < ct.min] <- 0
  
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }
  
  if (trace) cat('Begin GMPR size factor calculation ...\n')
  
  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm),  function(i) {		
    # if (i %% 50 == 0) {
    #   cat(i, '\n')
    # }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x / comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))		
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm=TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  }
  )
  
  if (sum(is.na(gmpr))) {
    warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'), 
                   '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                   'For these samples, their size factors are set to be NA! \n', 
                   'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                   'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
  }
  
  if (trace) cat('Completed!\n')
  if (trace) cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
  names(gmpr) <- names(comm.no) <- colnames(comm)
  
  attr(gmpr, 'NSS') <- comm.no
  
  return(gmpr)
}

mse <- function(x, y)
  mean((x - y) ^ 2)

analysis_wrapper <- function(Y, x, offset, 
                             method = c("edgeR", "DESeq2", "mgs")){
  
  q <- nrow(Y)
  S <- colSums(Y)
  if (method == "DESeq2"){
    
    xf <- factor(x)
    obj <-
      suppressMessages(
        DESeqDataSetFromMatrix(
          countData = Y,
          colData = data.frame(xf),
          design = ~ xf
        )
      )
    mynorms <- matrix(offset / S, byrow = T, nrow = nrow(Y), ncol = ncol(Y))
    DESeq2::normalizationFactors(obj) <- mynorms
    
    mod <- suppressMessages(results(DESeq(obj, quiet = TRUE)))
    beta1_hat <- log(2^mod@listData$log2FoldChange)
    pv <- mod@listData$pvalue
    
  } else if (method == "edgeR"){
    
    
    obj <- 
      edgeR::DGEList(
        counts = Y, 
        lib.size = offset,
        group = x
      )
    obj <- suppressMessages(edgeR::estimateDisp(obj))
    
    mod <- edgeR::exactTest(object = obj)
    
    beta1_hat <- log(2^mod$table$logFC)
    pv <- mod$table$PValue
    
    
  } else if (method == "PLN"){
    
    Yt <- t(Y)
    quiet(pln_out <- PLN(Yt ~ x + offset(log(offset)),
                         control = PLN_param(config_post =
                                               list(variational_var = F))))
    beta1_hat <- as.numeric(coef(pln_out)[2,])
    pv <- rep(NA, q)
    
    
  } else if (method == "mgs"){
    
    df <- data.frame(x = x)
    phenotypeData <-
      newMRexperiment(
        counts = Y,
        libSize = S,
        normFactors = offset,
      )
    
    mod <- fitFeatureModel(phenotypeData,
                           mod = model.matrix(~1 + x, df),
                           szero = T
    )
    out_dat <- MRcoefs(mod, number = q, group = 4)
    beta1_hat <- out_dat$logFC
    pv <- out_dat$pvalues
    
  }
  
  out <-
    data.frame(
      beta1_hat = beta1_hat,
      pv = pv
    )
  
  return(out)
  
}

gmean <- function(x){
  
  exp(mean(log(x)))
  
  
}

groupnorm <- function(Y, x, method = c("G-RLE", "FTSS"), 
                      prop_reference = 0.4,
                      use_median = FALSE
){
  
  method <- match.arg(method)
  if (!all(x %in% c(0,1))) stop ("Covariate must be binary.")
  n <- ncol(Y) # number of samples
  q <- nrow(Y) # number of taxa
  
  S <- colSums(Y)
  Yt <- t(Y)
  
  Y0 <- Yt[x == 0, ]
  Y1 <- Yt[x == 1, ]
  counts0 <- pmax(colSums(Y0), 1)
  counts1 <- pmax(colSums(Y1), 1)
  S0 <- sum(S[x == 0])
  S1 <- sum(S[x == 1])
  prop0 <- counts0 / S0
  prop1 <- counts1 / S1
  
  
  if (method == "G-RLE"){
    
    center <- exp(log(prop0)/2 + log(prop1)/2)
    scale0 <- median(prop0 / center)
    scale1 <- median(prop1 / center)
    norm0 <- scale0 * S0
    norm1 <- scale1 * S1
    offset0 <- (S / S0) * norm0
    offset1 <- (S / S1) * norm1
    offset <- ifelse(x==1, offset1, offset0)
    
  } else if (method == "FTSS"){
    
    log_ratios <- log(prop1 / prop0)
    if (use_median){
      lr_center <- median(log_ratios)
      
    } else{
      lr_dens <- density(log_ratios)
      lr_center <- with(lr_dens, x[which.max(y)])
    }
    
    center_quantile <- mean(log_ratios <= lr_center)
    q_lower <- max(center_quantile - prop_reference/2, 0)
    q_upper <- min(center_quantile + prop_reference/2, 1)
    lr_quants <- quantile(log_ratios, c(q_lower, q_upper))
    include <- as.numeric(log_ratios > lr_quants[1] &
                            log_ratios < lr_quants[2])
    include <- include * sum(S) / sum(include)
    offset <- as.numeric(Yt %*% include)
  }
  
  (offset / gmean(offset)) * gmean(S)
  
}

get_offset <- function(Y, x = NULL, method = c("TSS","CSS","RLE","TMM", "GMPR", "Wrench",
                                               "G-RLE", "FTSS")){
  method <- match.arg(method)
  Y_out <- Y
  n <- ncol(Y)
  p <- nrow(Y)
  libsize <- colSums(Y)
  
  if(is.null(rownames(Y))){
    rownames(Y) <- paste0("t", seq_len(p))
  }
  
  if (method == "TSS"){
    offset <- libsize
    
  } else if (method %in% c("TMM", "RLE")){
    Y <- apply(Y, 2, function(x) ifelse(x == 0, 1, x)) # pseudo-count
    sf <- edgeR::calcNormFactors(Y, method = method) # scaling factor
    offset <- sf * libsize
    
  } else if (method == "CSS"){
    
    Y_exp <- newMRexperiment(Y)
    suppressMessages(perc <- cumNormStatFast(Y_exp))
    sf <- metagenomeSeq::calcNormFactors(Y_exp, p = perc)
    sf <- sf$normFactors / libsize 
    sf <- sf / exp(mean(log(sf)))
    offset <- sf * libsize
    
  } else if (method == "GMPR"){
    
    sf <- as.numeric(gmpr(Y, trace = F)) / libsize
    sf <- sf / gmean(sf)
    offset <- sf * libsize
    
  } else if (method == "Wrench"){
    
    xf <- as.factor(x)
    w_out <- Wrench::wrench(Y, condition = xf)
    offset <- as.vector(w_out$nf)
    offset <- offset * gmean(libsize) / gmean(offset)
    
  } else {
    
    if (is.null(x)) stop ("Must provide covariate to use group-wise normalization.")
    offset <- groupnorm(Y = Y, x = x, method = method)  
    
  }
  
  offset
  
}


get_synthetic_data <- function(Y, beta1, S){
  
  
  # Set Up
  n <- ncol(Y)
  n0 <- ceiling(n/2)
  n1 <- n - n0
  x <- c(rep(0,n0),rep(1,n1))
  x1_id <- n0 + seq_len(n1)
  fold_change <- exp(beta1)
  
  # Shuffle Data
  perm <- sample(n)
  Y <- Y[, perm]
  S <- S[perm]
  
  # Re-sampling
  for (i in x1_id) {
    Y[, i] <-
      rmultinom(
        n = 1,
        size = S[i],
        prob = (Y[, i] / S[i]) * fold_change
      )
  }

  list(Y = Y, x = x, S = S)
}

# Simulation functions ----------------------------------------------------

general_model_sims <-
  function(settings, reps_do = 1, seed_param = 1234){
    
    # Settings
    nsettings <- nrow(settings)
    n_list <- settings$n # sample size
    s2_v_list <- settings$s2_v # variance
    ps_list <- settings$ps # prop signal
    norm_list <- settings$norm # norm method
    method_list <- settings$method # DAA method
    
    if ("signal" %in% colnames(settings)){
      
      signal_list <- settings$signal
      
    } else {
      
      signal_list <- rep("right", nsettings)
      
    }
    
    
    if ("libscale" %in% colnames(settings)){
      libscale_list <- settings$libscale # group1 library size scalar
      
    } else {
      libscale_list <- rep(1, nsettings)
      
    }
    
    
    # Non-varied parameters
    a0 <- -qnorm(c(0.5))
    q <- 300
    mean_lib <- 60000
    var_lib <- 15000^2
    size_lib <- mean_lib^2 / (var_lib - mean_lib)
    
    fdr_results <- list()
    fpr_results <- list()
    j <- 1
    
    for (s in seq_len(nsettings)){
      
      # Settings
      ps <- ps_list[s]
      n <- n_list[s]
      s2_v <- s2_v_list[s]
      norm <- norm_list[s]
      method <- method_list[s]
      signal <- signal_list[s]
      libscale <- libscale_list[s]
      
      # Sample signals
      set.seed(seed_param)
      beta0 <- rnorm(q, sd = 1)
      q1 <- ceiling(q * ps)
      beta1 <- rep(0, q)
      true_beta <- logical(q)
      if (q1 > 0){
        
        temp <- rnorm(q1)
        temp <- temp - mean(temp)
        if (signal == "right"){
          
          beta1[seq_len(q1)] <- temp + 1
          
        } else if (signal == "center"){
          
          beta1[seq_len(q1)] <- temp
          
        } else if (signal == "mixed"){
          
          beta1[seq_len(q1)] <- temp + ifelse(seq_len(q1) <= q1/2, -1, 1)
          
        }
        true_beta[seq_len(q1)] <- TRUE
        
      }
      
      # Sample covariance matrices
      R <-
        cov2cor(LaplacesDemon::rinvwishart(2*q, diag(q)))
      Sigma <- R * s2_v
      
      # Generate X, nu
      x <- numeric(n) + 1
      x[seq_len(n/2)] <- 0
      nu <- exp(beta0 + beta1 %*% t(x))
      
      for (r in reps_do){
        
        # Generate data
        set.seed(r + 1000)
        S <- (rnbinom(n = n, size = size_lib, mu = mean_lib))
        S[x==1] <- floor(S[x==1] * libscale)
        
        # Generate latent variables
        W <- t(mvtnorm::rmvnorm(n = n, sigma = R))
        V <- t(mvtnorm::rmvnorm(n = n, sigma = Sigma))
        not_zero <- W + a0 > 0
        lambda <- nu * exp(V)
        mu <- lambda * not_zero
        
        # Sample data
        Y <- matrix(NA, q, n)
        for (i in 1:n) {
          Y[, i] <- as.numeric(rmultinom(1, S[i], prob = mu[,i]))
          
        }
        
        Yt <- t(Y)
        
        
        if (norm != "linda"){
          # Calculate normalization factor
          quiet(offset <- get_offset(Y = Y, x = x, method = norm))
          
          # Analyze
          out <- analysis_wrapper(Y, x, offset, method = method)
          beta1_hat <- out$beta1_hat
          pv <- out$pv
          
        } else {
          
          mod <- 
            quiet(MicrobiomeStat::linda(
              feature.dat = Y,
              meta.dat = data.frame(x = x),
              formula = "~ x",
              feature.dat.type = "count",
              # zero.handling = "pseudo-count",
              adaptive = T
            ))
          
          out <-
            with(mod$output$x,
                 data.frame(beta1_hat = log(2 ^ log2FoldChange),
                            pv = pvalue)
            )
        }
        
        
        # Calculate MSE
        mse_all <- mse(beta1, beta1_hat)
        mse0 <- mse(beta1[!true_beta], beta1_hat[!true_beta])
        mse1 <- NA
        if (q1 > 0){
          mse1 <- mse(beta1[true_beta], beta1_hat[true_beta])
        }
        
        # FDR results 
        cutoffs <- c(0.05, 0.10)
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
            signal = signal,
            libscale = libscale,
            s2_v = s2_v,
            q = q,
            q1 = q1,
            fdr = cutoffs, # true FDR
            tpr = tpr, 
            fpr = fpr,
            fdp = fdp,
            mse_all = mse_all,
            mse0 = mse0,
            mse1 = mse1
          )
        
        
        # FPR/TPR results
        false_pv_indices <- which(!(true_beta[order(pv)]))
        fp <- seq_len(q - q1)
        fpr <- fp / (q- q1)
        tp <- false_pv_indices - fp
        tpr <- tp / q1
        
        df <-
          data.frame(
          method = method,
          norm = norm,
          n = n,
          r = r,
          signal = signal,
          libscale = libscale,
          s2_v = s2_v,
          q = q,
          q1 = q1,
          tpr = tpr, 
          fpr = fpr,
          tp = tp,
          fp = fp
        )
        
        fpr_results[[j]] <- subset(df, df$fpr <= 0.40)

        j <- j + 1
        
      }
      
    }
    
    fdr_results_sum <- bind_rows(fdr_results)
    fpr_results_sum <- bind_rows(fpr_results)
    list(fdr = fdr_results_sum,
         fpr = fpr_results_sum)
    
    
    
  }

general_data_sims <-
  function(settings, reps_do = 1, seed_param = 1234){
    
    # Load datasets
    dataset_lists <- list()
    load("data/PHACS.rda")
    dataset_lists$PHACS <- 
      list(
        Y = phacs$taxa,
        S = colSums(phacs$taxa)
      )
    load("data/MLVS_MBS.rda")
    dataset_lists$HPFS <-
      list(
        Y = mlvs_mbs$taxa,
        S = colSums(mlvs_mbs$taxa)
      )
    
    # Settings
    nsettings <- nrow(settings)
    data_list <- settings$data # sample size
    ps_list <- settings$ps # prop signal
    norm_list <- settings$norm # norm method
    method_list <- settings$method # DAA method
    
    if ("signal" %in% colnames(settings)){
      
      signal_list <- settings$signal
      
    } else {
      
      signal_list <- rep("right", nsettings)
      
    }
    
    fdr_results <- list()
    fpr_results <- list()
    j <- 1
    for (s in seq_len(nsettings)){
      
      # Settings
      ps <- ps_list[s]
      norm <- norm_list[s]
      method <- method_list[s]
      signal <- signal_list[s]
      data <- data_list[s]
      
      # Load data
      set.seed(1)
      Y_true <- dataset_lists[[data]]$Y
      S <- dataset_lists[[data]]$S
      q <- nrow(Y_true)
      n <- ncol(Y_true)
      Y_true <- Y_true[sample(q),]
      
      # Sample signals
      set.seed(seed_param)
      beta0 <- rnorm(q, sd = 1)
      q1 <- ceiling(q * ps)
      beta1 <- rep(0, q)
      true_beta <- logical(q)
      if (q1 > 0){
        
        temp <- rnorm(q1)
        temp <- temp - mean(temp)
        if (signal == "right"){
          
          beta1[seq_len(q1)] <- temp + 1
          
        } else if (signal == "center"){
          
          beta1[seq_len(q1)] <- temp
          
        } else if (signal == "mixed"){
          
          beta1[seq_len(q1)] <- temp + ifelse(seq_len(q1) <= q1/2, -1, 1)
          
        }
        true_beta[seq_len(q1)] <- TRUE
        
      }
      
      for (r in reps_do){
        
        # Generate data
        set.seed(r + 1000)
        synth_data <- get_synthetic_data(Y = Y_true, beta1 = beta1, S = S)
        Y <- synth_data$Y
        x <- synth_data$x
        Yt <- t(Y)
        
        if (norm != "linda"){
          # Calculate normalization factor
          quiet(offset <- get_offset(Y = Y, x = x, method = norm))
          
          # Analyze
          out <- analysis_wrapper(Y, x, offset, method = method)
          beta1_hat <- out$beta1_hat
          pv <- out$pv
          pv[is.na(pv)] <- 1
        } else {
          
          mod <- 
            quiet(MicrobiomeStat::linda(
              feature.dat = Y,
              meta.dat = data.frame(x = x),
              formula = "~ x",
              feature.dat.type = "count",
              # zero.handling = "pseudo-count",
              adaptive = T
            ))
          
          out <-
            with(mod$output$x,
                 data.frame(beta1_hat = log(2 ^ log2FoldChange),
                            pv = pvalue)
            )
        }
        
        
        # Calculate MSE
        mse_all <- mse(beta1, beta1_hat)
        mse0 <- mse(beta1[!true_beta], beta1_hat[!true_beta])
        mse1 <- NA
        if (q1 > 0){
          mse1 <- mse(beta1[true_beta], beta1_hat[true_beta])
        }
        
        # FDR results 
        cutoffs <- c(0.05, 0.10)
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
            data = data,
            r = r,
            signal = signal,
            q = q,
            q1 = q1,
            fdr = cutoffs, # true FDR
            tpr = tpr, 
            fpr = fpr,
            fdp = fdp,
            mse_all = mse_all,
            mse0 = mse0,
            mse1 = mse1
          )
        
        
        # FPR/TPR results
        false_pv_indices <- which(!(true_beta[order(pv)]))
        fp <- seq_len(q - q1)
        fpr <- fp / (q- q1)
        tp <- false_pv_indices - fp
        tpr <- tp / q1
        

        
        df <-
          data.frame(
            method = method,
            norm = norm,
            data = data,
            r = r,
            signal = signal,
            q = q,
            q1 = q1,
            q = q,
            q1 = q1,
            tpr = tpr, 
            fpr = fpr,
            tp = tp,
            fp = fp
          )
        
        fpr_results[[j]] <- subset(df, df$fpr <= 0.40)
        j <- j + 1
        
      }
      
    }
    
    fdr_results_sum <- bind_rows(fdr_results)
    fpr_results_sum <- bind_rows(fpr_results)
    list(fdr = fdr_results_sum,
         fpr = fpr_results_sum)
    
    
    
  }




