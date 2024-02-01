########################################
# XCMS - PART 1: Peak identification 
########################################
# https://bioconductor.org/packages/devel/bioc/vignettes/xcms/inst/doc/xcms.html
source("~/Git/exAtlas/src/R/XCMS_pipeline/XCMS_functions.R")
suppressPackageStartupMessages({
  library(magrittr)
})
tictoc::tic()

################################################           
# Load files  
################################################
# Output directory
out_dir <- "~/Sara/exAtlas/res/XCMS/MS"
in_file <- "~/Sara/exAtlas/meta/exAtlas_MS_meta_info.tsv"
  
# Output name 
file_name <- "exAtlas_MS" 

# Meta data file
pd <- read.delim(in_file)

################################################           
# Set XCMS parameters
################################################
# Data frame with meta data
pd = pd 

# Centwave Parameters (use AutoTuner to help set these)
peakwidth = c(2, 20) 
ppm = 10
noise = 120
snthresh = 3
integrate = 2
prefilter = c(2,250)
mzdiff = -0.001
fitgauss = FALSE 
mzCenterFun = "wMean"

# Filtering: Intensity Parameter 
ReFine = TRUE # if TRUE, data will be filtered
thr = 3000 # set intensity threshold (i.e., everything <X will be removed)

# Alignment: Peak Groups Parameters 
Align = FALSE # if TRUE, data will be aligned
minFraction = 0.85
span = 0.4
smooth = "linear"
subsetAdjust = "average"
sampleGroupsType = pd$Type # definition of the sample groups (e.g., samples, QC, Blank).
sampleGroupsType_sel = "QC" # the subset of samples on which the alignment of the full data set will be based (e.g. subset being the index of QC samples)

################################################           
# Check assumptions
################################################
# Check if downstream dataframe columns exists  
if(!"File" %in% colnames(pd) | !"Genotype" %in% colnames(pd) | !"Genotype_nr" %in% colnames(pd) | !"Injection" %in% colnames(pd) | !"Sample" %in% colnames(pd)){
  stop("Column(s) missing in pd, fix column(s) or adapt script")
}

# Check if output directory exists 
if(!dir.exists(out_dir)){
  stop(paste("Output directory does not exist:", out_dir))
} else{
  message(paste("Input file:", in_file))
  dir.create(file.path(out_dir, "Peak_identification"), showWarnings = FALSE)
  message(paste("Output directory:", file.path(out_dir, "Peak_identification")))
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
res_PeakID <- XCMS_peak_identification(pd, 
                                       peakwidth = peakwidth, 
                                       ppm = ppm, 
                                       noise = noise, 
                                       snthresh = snthresh, 
                                       integrate = integrate, 
                                       prefilter = prefilter, 
                                       mzdiff = mzdiff, 
                                       fitgauss = fitgauss,
                                       mzCenterFun = mzCenterFun,
                                       ReFine = ReFine,
                                       thr = thr, 
                                       sampleGroupsType = sampleGroupsType,
                                       sampleGroupsType_sel = sampleGroupsType_sel,
                                       Align = Align,
                                       minFraction = minFraction,
                                       span = span,
                                       subsetAdjust = subsetAdjust,
                                       smooth = smooth,
                                       nCore = nCore)

# Save data 
save(res_PeakID, file = file.path(out_dir, "Peak_identification", paste0(file_name, "_XCMS_peak_identification.RData")))

# Save session info
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "Peak_identification", paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")
tictoc::toc()