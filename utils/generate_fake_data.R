
current_objects <- ls()
set.seed(1)

# PHACS
q <- 344
n <- 254

beta0 <- rnorm(q)
x <- rep(c(0,1), each = n/2)
beta1 <- rep(0,q)
q1 <- floor(q * 0.1)
beta1[seq_len(q1)] <- rnorm(q1, mean = 1)
s2_v <- 1.5
Sigma_z <- cov2cor(LaplacesDemon::rinvwishart(2 * q, diag(q)))
Sigma_v <- Sigma_z * s2_v
Z <- t(mvtnorm::rmvnorm(n = n, sigma = Sigma_z))
V <- t(mvtnorm::rmvnorm(n = n, sigma = Sigma_v))
mu <- exp(beta0 + V) * (Z > 0) * exp(beta1 %*% t(x))
mean_lib <- 60000
var_lib <- 15000^2
size_lib <- mean_lib^2 / (var_lib - mean_lib)
S <- rnbinom(n = n, size = size_lib, mu = mean_lib)
Y <- matrix(NA, nrow = q, ncol = n)
for (i in seq_len(n)) Y[, i] <- as.numeric(rmultinom(1, S[i], prob = mu[, i]))

rownames(Y) <- paste0("t",seq_len(q))
phacs <- list()
phacs$taxa <- Y
phacs$x <- x
save(phacs, file = "data/PHACS.rda")

# MLVS/MBS
q <- 372
n <- 520

beta0 <- rnorm(q)
x <- rep(c(0,1), each = n/2)
beta1 <- rep(0,q)
q1 <- floor(q * 0.1)
beta1[seq_len(q1)] <- rnorm(q1, mean = 1)
s2_v <- 1.5
Sigma_z <- cov2cor(LaplacesDemon::rinvwishart(2 * q, diag(q)))
Sigma_v <- Sigma_z * s2_v
Z <- t(mvtnorm::rmvnorm(n = n, sigma = Sigma_z))
V <- t(mvtnorm::rmvnorm(n = n, sigma = Sigma_v))
mu <- exp(beta0 + V) * (Z > 0) * exp(beta1 %*% t(x))
mean_lib <- 60000
var_lib <- 15000^2
size_lib <- mean_lib^2 / (var_lib - mean_lib)
S <- rnbinom(n = n, size = size_lib, mu = mean_lib)
Y <- matrix(NA, nrow = q, ncol = n)
for (i in seq_len(n)) Y[, i] <- as.numeric(rmultinom(1, S[i], prob = mu[, i]))
rownames(Y) <- paste0("t",seq_len(q))
mlvs_mbs <- list()
mlvs_mbs$taxa <- Y
mlvs_mbs$x <- x
save(mlvs_mbs, file = "data/MLVS_MBS.rda")

rm(list = setdiff(ls(), current_objects))


