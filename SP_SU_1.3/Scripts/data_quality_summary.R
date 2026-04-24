#####################################################################
#====================================================================
# Script Name : data_type_plots_WSP2024_chinook.R
# Script Purpose: Create estimate class plots for FIA Chinook Spring and
#                 Summer 1.3 CUs
#                 to include in the WSP rapid status templates. 
# Script Date : March 30 2026
# Script Author: Isabella Borea 
# R Version : R version 4.2.2 (2024-10-31 ucrt)
# R Packages : tidyverse, readxl, ggplot2
# Control + Shift + C = Comment selected lines
#====================================================================
#####################################################################

library(tidyverse)
library(readxl)
library(ggplot2)
library(lubridate)
library(scales)
library(RColorBrewer)
library(dplyr)

###########################################

quality_data <- read.csv("NuSEDS_CUdata_1.3_FORMATTED_08212025_OUT", strip.white = TRUE)

str(quality_data)
glimpse(quality_data)
as.tibble(quality_data)





