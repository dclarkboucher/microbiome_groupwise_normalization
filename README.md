# microbiome_groupwise_normalization
 A code and data repository for the article "Group-wise normalization in differential abundance analysis of microbiome samples"
 
 The first step in using this repository is to open the file "microbiome_groupwise_normalization.Rproj", an R project that 
 will automatically open RStudio in the correct working directory. The full analysis comprises three categories of simulations: 
 model-based simulations based on the multinomial distribution, synthetic data simulations based on real microbiome count 
 matrices, and model-based simulations with a confounding variable. Once the R project has been opened, these analyses can be 
 implemented with the individuals scripts "model_based_simulations.R", "synthetic_data_simulations.R", and 
 "model_based_simulations_confounder.R", respectively, which are located in the "analysis" directory. In theory, the whole 
 analysis could be replicated by setting the "total_reps" parameter in these scripts to 1,000, telling R to run a thousand 
 replications of each setting. However, this selection would be intensely computationally heavy and would likely require 
 parallelization on a remote computing cluster. Instead, to replicate a fraction of the study, one could set the "total_reps" 
 parameters to say, 50, which may be more manageable. In addition, while the real microbiome datasets used in the synthetic data 
 analysis are not currently included in this repository, the "synthetic_data_simulations.R" script will automatically generate a 
 toy dataset that can be analyzed for demonstrative purposes. After running the analysis locally, the figures can be generated 
 with the script "analysis/process_output.R", and will be saved in the "figures" directory. 
