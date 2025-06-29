install.packages(c("tibble", "dplyr", "MicrobiomeStat", "LaplacesDemon",
                   "mvtnorm", "tidyr", "forcats", "ggpubr"))
install.packages("BiocManager") # need this for installing other packages

BiocManager::install(c("edgeR", "DESeq2","metagenomeSeq"))
