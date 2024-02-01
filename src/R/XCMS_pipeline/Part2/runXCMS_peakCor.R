#######################################
# XCMS - Part 2: Peak correspondance
#######################################
# https://gitlab.com/stanstrup-teaching/XCMS-course/-/blob/master/1-XCMS.Rmd
source("~/Git/exAtlas/src/R/XCMS_pipeline/XCMS_functions.R")
suppressPackageStartupMessages({
  library(xcms)
  library(magrittr)
  })
tictoc::tic()

################################################           
# Load files  
################################################
# Output directory 
out_dir <- "~/Sara/exAtlas/res/XCMS/MS"
in_file <- "~/Sara/exAtlas/res/XCMS/MS/Peak_identification/exAtlas_MS_XCMS_peak_identification.RData"
  
# Output file prefix
file_name <- "exAtlas_MS"

# XCMS peakID data
load(in_file)

################################################           
# Set XCMS parameters. Use AutoTuner for support. 
################################################
# Set xcms object
XCMS_obj = "Peak_data_refined" 
peakID_obj = res_PeakID[[XCMS_obj]]

# Set PeakDensity parameters. Tips! To help set these, look at the link at the top of this page 
sampleGroups = pData(peakID_obj)$Tissue_exp  
minFraction = 0.6  
bw = 3.276
binSize = 0.01
minSamples = 2
maxFeatures = 50

# Set quantify parameter
value = "maxo"

# Should the peaks be filled? 
PeakFill = TRUE

################################################           
# Free some space + sanity check
################################################
message(paste("XCMS object selected for correspondance:", XCMS_obj))
rm(res_PeakID)

# Check if output directory exist 
if(!dir.exists(out_dir)){
  stop(paste("Output directory does not exist:", out_dir))
} else{
  message(paste("Input file:", in_file))
  dir.create(file.path(out_dir, "Peak_correspondance"), showWarnings = FALSE)
  message(paste("Output directory:", file.path(out_dir, "Peak_correspondance")))
  message(paste("Output file prefix:", file_name))
}

# Set parallelisation - nr of cores
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1) {
  stop("Nr of cores for parallelisation must be supplied", call.=FALSE)
} else {
  nCore <- args[1]
  message("Nr of cores: ", nCore)
}

################################################           
# Run XCMS 
################################################
res_cor <- XCMS_correspondance(peakID_obj = peakID_obj, 
                               sampleGroups = sampleGroups, 
                               minFraction = minFraction, 
                               bw = bw, 
                               binSize = binSize, 
                               minSamples = minSamples, 
                               maxFeatures = maxFeatures, 
                               value = value,
                               PeakFill = PeakFill, 
                               nCore = nCore)

################################################           
# Save data 
################################################
save(res_cor, file = file.path(out_dir, "Peak_correspondance", paste0(file_name, "_XCMS_peak_correspondance.RData")))

# Save session info
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "Peak_correspondance", paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")
tictoc::toc()