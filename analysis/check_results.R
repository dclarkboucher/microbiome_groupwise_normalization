
load("output/model_sims.rda")
fdr <- results_sum$fdr
dim(fdr)

fdr1 <- 
  fdr %>% 
  filter(
    r == 1,
    round(fdr,2) == 0.10,
    norm == "TSS",
    method == "metagenomeSeq",
    q1 == 90,
    n == 200
  ); fdr1
