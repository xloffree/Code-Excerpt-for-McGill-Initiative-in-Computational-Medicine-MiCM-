# Code Excerpt for McGill Initiative in Computational Medicine (MiCM) 
Conduct a Wilcoxon test on the before and after effector index (EI) scores of the genes associated with each single nucleotide polymorphism (SNP). The before EI scores are the originals. The after EI scores are calculated by running the algorithm iteratively with each SNP removed. The Wilcoxon test is meant to determine how much the order of the EI scores associated with the genes associated with each SNP changes.

To run, simply highlight the entire script titled "micm2.R" and click run. The csv file titled "original_and_new_EI_calcium2.csv" is input for the script and is required to run the script. The desired results for the script can be found in the csv file titled "wcx_snp_calcium.csv".
