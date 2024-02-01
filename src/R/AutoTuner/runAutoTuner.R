suppressPackageStartupMessages({
  library(Autotuner)
  library(dplyr)
  library(tidyverse)
})
tictoc::tic()

################################################           
# Load files  
################################################
pd <- read_tsv("~/Sara/exAtlas/meta/exAtlas_MS_meta_info.tsv", show_col_types = FALSE)
file_name <- "exAtlas_subsamples_Autotuner_res"
out_dir <- "~/Sara/exAtlas/res/AutoTuner"

# OBS! I recommend that you use SLURM for 'Parameter Extraction from Individual Extracted Ion Chromatograms'. The rest should, if possible, be interactive.
################################################           
# Check assumptions before proceeding 
################################################
# Check if output directory exists 
if(!dir.exists(out_dir)){
  stop(paste("Output directory does not exist:", out_dir))
} else{
  message(paste("Output directory:", out_dir))
  message(paste("Save output file containing:", file_name))
}

########################################################           
# Set parameters for AutoTuner Sliding Window Analysis
########################################################
# The file_col argument corresponds to the string column of the metadata that matches raw data samples by name. 
# The factorCol argument corresponds to the metadata column containing information about which experimental factor each sample belongs to. 
Autotuner <- createAutotuner(pd$File,
                             pd,
                             file_col = "Sample",
                             factorCol = "Tissue")

# The user should prioritize finding where peaks *start* rather than caturing the entire peak bound. 
# Downstream steps actually do a better job of estimating what the proper peak bounds should be. 
# The user should play with the lag, threshold, and influence parameters to perform the sliding window analysis. 
# Lag - The number of chromatographic scan points used to test if next point is significance (ie the size number of points making up the moving average). 
# Threshold - A numerical constant representing how many times greater the intensity of an adjacent scan has to be from the scans in the sliding window to be considered significant. 
# Influence - A numerical factor used to scale the magnitude of a significant scan once it has been added to the sliding window. 
lag <- 20
threshold<- 1
influence <- 0.3

########################################################           
# Total Ion Current Peak Identification 
########################################################
message("Sliding window analysis")
AutoTuner_param_testing <- list(Lag = lag,
                                Threshold = threshold,
                                Influence = influence) 

signals <- lapply(getAutoIntensity(Autotuner), 
                  ThresholdingAlgo, lag, threshold, influence)

pdf(file.path(out_dir, paste0(file_name, "_sliding_window_res.pdf")))
plot_signals(Autotuner, 
             threshold, 
             ## index for which data files should be displayed
             sample_index = 1:nrow(pd), 
             signals = signals)
graphics.off()
rm(lag, influence, threshold)

### Interpreting Sliding Window Results ###
# The figure above has two components:
# 1) Top Plot: The chromatotgraphic trace for each sample (solid line) along with the noise associated with each sample (dashed line). 
# 2) Bottom Plot: A signal plot used to indicate which chromatographic regions have peaks.
# The user should look for combinations of the three sliding window parameters that returns many narrow peaks within the signal plot. See the example above.
# Autotuner will expand each of these regions to obtain improved estimates on the bounds within the isolatePeaks function below.

########################################################           
# Checking Peak Estimates
########################################################
# The peaks with expanded bounds returned from the isolatePeaks function can be rapidly checked visually using the plot_peaks function as shown below. 
# The bounds should capture the correct ascention and descention points of each peak. If peak bounds are not satisfactory, the user should return to the sliding window analysis, and try a different conbination of the three parameters.
# Remember, this whole process is only designed to isolate regions enriched in real features rather than find true peaks. The bounds don't need to be completely perfect. 
# Its much more important that the bounds contain some kind of chromatographic peaks rather than less dynamic regions of the chromatographic trace.
Autotuner <- isolatePeaks(Autotuner = Autotuner, 
                          returned_peaks = 10, 
                          signals = signals)

pdf(file.path(out_dir, paste0(file_name,"_peak_estimate_check.pdf")))
for(i in 1:50) {
  plot_peaks(Autotuner = Autotuner, 
             boundary = 150, 
             peak = i)    
}
graphics.off()

#####################################################################           
# Parameter Extraction from Individual Extracted Ion Chromatograms
#####################################################################
# In order to estimate parameters from the raw data, the user should run the EICparams function as below. 
# The massThreshold is an absolute mass error that should be greater than the expected analytical capabilities of the mass analyzer.
message("Parameter extracion from EIC")
eicParamEsts <- EICparams(Autotuner = Autotuner, 
                          massThresh = .005, 
                          verbose = FALSE,
                          returnPpmPlots = FALSE,
                          useGap = TRUE)

#####################################################################           
# Save data 
#####################################################################
# All that remains now is to get what the dataset estimates are.
# Running AutoTuner is now complete, and the estimates may be entered directly into XCMS to processes raw untargeted metabolomics data
message("Save data")
xcms_param <- returnParams(eicParamEsts, Autotuner)
save(eicParamEsts, Autotuner, xcms_param, AutoTuner_param_testing, file = file.path(out_dir, paste0(file_name, ".RData")))

# Session Info
# Save session info
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")
tictoc::toc()