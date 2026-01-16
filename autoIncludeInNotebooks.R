
#' Packages and source code to autoinclude in all notebooks
cat("Load Packages...")
# Favorite data wrangling and functionality packages
suppressPackageStartupMessages(expr = {
library(data.table) 
library(DT)
library(pbapply)
library(doParallel)
library(readxl)
library(rstudioapi)
library(dynamicTreeCut)

# Favorite data vis packages
library(circlize)
library(viridis)
library(ggplot2)
library(ggrepel)
library(eulerr)
library(ComplexHeatmap)
})

# Favorite Source Code
source("~/Documents/Code/bp_utils/ManageScriptData.R")
source("~/Documents/Code/bp_utils/UniprotIDMapping.R")

source("~/Documents/Code/maxb_utils/GenericHelpers.R")

#' Names and file locations of common optional packages/scripts:

#' source("~/Documents/Code/maxb_utils/BioinformaticsHelpers.R")
#' source("~/Documents/Code/maxb_utils/ProteomicsPlots.R")
#' source("~/Documents/Code/bp_utils/SaintExpressR.R")
#' source("~/Documents/Code/bp_utils/MSstats_V4_Functions.R")
#' source("~/Documents/Code/bp_utils/MSstats_Helper_Functions.R")
#' source("~/Documents/Code/bp_utils/enrichmentTestFunctions.R")
#' source("~/Documents/Code/bp_utils/NMF_Helpers.R")
#' source("~/Documents/Code/bp_utils/APMS_Utils.R")
         
#'   Bioinformatic Workflow Packages
#' library(artMS)
#' library(MSstats)
#' library(MSstatsConvert)
#' library(cRomppass)
#' library(clusterProfiler)
#' library(fgsea)
         
if (!file.exists(paste0(getwd(),"/logs"))){
  dir.create(paste0(getwd(),"/logs"))
}
dataFolder <- ScriptNamedDir()
calledFrom <- substrLeft( tail(strsplit(dataFolder,"/") [[1]], n = 1), 0, -5)

cat("Log installed package versions...")
WriteSessionInfo(path = paste0( "./logs/", DateFileName(paste0(calledFrom, "_SessionInfo.txt")) ) )
WriteInstalledPackages(path = paste0( "./logs/", DateFileName(paste0(calledFrom, "_InstalledPackages.txt")) ) )

