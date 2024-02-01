suppressPackageStartupMessages({
  library(RAMClustR)
  library(readr)
  library(tidyverse)
})
tictoc::tic()

################################################           
# Load data
################################################
out_dir <- "~/Sara/exAtlas/res/RAMClustR/MS"
in_file <- "~/Sara/exAtlas/res/XCMS/MS/Peak_correspondance/exAtlas_MS_XCMS_peak_correspondance.RData"
file_name <- "exAtlas_MS"

load(in_file)

################################################           
# Set RAMClustR parameters 
################################################
XCMS_obj <- res_cor$PeakID_obj_gr_filled
normalize_type = "none"
minCluster = 1 
sigma_t = NULL
sigma_r = NULL
max_timedif = NULL

mode = "negative"
mzabs.error = 0.02
ppm.error = 10

################################################           
# Check assumptions before proceeding 
################################################
if(!dir.exists(out_dir)){
  stop(paste("Output directory does not exist:", out_dir))
} else{
  message(paste("Input file:", in_file))
  setwd(out_dir)
  message(paste("Output directory:", out_dir))
  message(paste("Output file prefix:", file_name))
}

################################################           
# Run RAMClustR
################################################
# First column in meta data (see pData(XCMS_obj)) will be used as sample labels 
message("Run RAMClustR")
experiment <- defineExperiment(force.skip = TRUE)
RC_res <- ramclustR(xcmsObj = XCMS_obj, 
                    ExpDes = experiment,
                    normalize = normalize_type,
                    minModuleSize = minCluster, 
                    st = sigma_t,
                    sr = sigma_r,
                    maxt = max_timedif
                    )

################################################           
# Add ion info 
################################################
message("Run interpret MS spectra function")
RC_res <- do.findmain(RC_res, 
                      mode = "negative", 
                      mzabs.error = mzabs.error, 
                      ppm.error = ppm.error
                      )

# Save RC settings
RC_param <- list(
  "normalize" = normalize_type,
  "mode" = mode,
  "mzabs_error" = mzabs.error,
  "ppm_error" = ppm.error,
  "minModuleSize" = minCluster,
  "st" = sigma_t,
  "sr" = sigma_r,
  "maxt" = max_timedif
  )

################################################           
# Save data 
################################################
message(paste("Saving data"))
save(RC_res, RC_param, file = paste0(file_name, "_RAMClust_res.RData"))

writeLines(capture.output(sessionInfo()), con = paste("sessionInfo",Sys.Date() ,".txt", sep=""))
message("Done")

tictoc::toc()