### FIA-CN-FSRR-FSAR ###
General Fraser and Interior Area Chinook repository for FSRR/FSAR analyses



# Instructions for reproducing the Fraser Spring 1.2 Chinook Analysis
The bulk of the scripts and data for reproducing this analysis are housed in the "SP_1.2" folder. The folder is organized into "Data_in", "Figures", "Results", and "Scripts". Most data referenced for the analysis can be found in "Data_in". Occasionally, data we anticipate will be useful for future analyses for other SMUs (Fraser Spring and Summer 1.3 Chinook) are kept in the folder "data" within the main repository. Similarly, scripts that are likely to be re-used for future SMU analyses are kept in the "R" folder within the main repository.

To re-run the analysis, go into the "SP_1.2" folder and open the "Scripts" sub-folder. You will find a series of R-Markdown files that are numbered in the order they need to be run to successfully reproduce the analysis. As of 27-06-2025, this series of scripts function in the order they are laid out using R V4.4.2 and Stan V2.32.2.Please note that there is a specific version of the "samSim" package that has to be installed directly from Github - the details are provided in the first .Rmd where this package is required. Finally, the results of the analysis will vary slightly from run to run because of the Bayesian approach to the analysis. 