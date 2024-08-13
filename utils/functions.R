library(matrixStats)
library(metagenomeSeq)
# packages: edgeR, 

# Function to remove unwanted output
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Function to perform DAA
# Y: q by n matrix of q taxa, n samples
# x: Binary covariate indicating group status, length n
# offset: offset term of length n
# method: DAA method (edgeR, DESeq2, metagenomeSeq)
analysis_wrapper <- function(Y, x, offset, 
                             method = c("edgeR", "DESeq2", "metagenomeSeq")){
  
  q <- nrow(Y)
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
    mynorms <- matrix(offset / colSums(Y), byrow = T, nrow = nrow(Y), ncol = ncol(Y))
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
    
    
  } else if (method == "metagenomeSeq"){
    
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





# Function to compute the geometric mean
gmean <- function(x) exp(mean(log(x)))


# Function to perform group-wise normalization
# Y: q by n matrix of q taxa, n samples
# x: Binary covariate indicating group status
# method: the normalization method, G-RLE or FTSS
# prop_reference: proportion of taxa to be used as reference taxa in FTSS
# use_median: flag to use median in FTSS instead of mode
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

# Function to calculate the offset term
# Y: q by n matrix of q taxa, n samples
# x: (optional) Binary covariate indicating group status, length n
# method: normalization method
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

gmean <- function(x){
  
  exp(mean(log(x)))
  
  
}

# Reference: https://github.com/lichen-lab/GMPR
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





