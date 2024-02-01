suppressPackageStartupMessages({
  library(tidyverse)
  library(xcms)
})

# Load data
load("~/Sara/exAtlas/res/XCMS/MS/Peak_correspondance/exAtlas_MS_XCMS_peak_correspondance.RData")
load("~/Sara/exAtlas/res/RAMClustR/MS/exAtlas_MS_RAMClust_res.RData")

# Set output 
out_dir <- "~/Sara/exAtlas/res/Preprocessing_res/MS"
file_name <- "exAtlas_MS"

RC_info <- 
  data.frame(mz = RC_res$M, 
             rt = RC_res$clrt) %>% 
  rownames_to_column("RC")

# Prep RAMClustR table by extracting columns of interest from RC_res
RC_cluster_res <- RC_res$SpecAbund %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(mz = RC_res$M, 
         rt = RC_res$clrt, 
         pre_ion = RC_res$precursor.mz, 
         pre_ion_type = RC_res$precursor.type) %>% 
  rownames_to_column("RC") %>% 
  mutate(RC = str_remove(str_remove(RC, "^C+"), "^0+") %>% as.integer()) %>% 
  select(c(RC, mz, rt, pre_ion, pre_ion_type, everything()))

# Extract xcms peak info (full and refined version)
xcms_peak_info <- featureDefinitions(res_cor$PeakID_obj_gr_filled) %>% 
  as.data.frame() %>% 
  rownames_to_column("Feature") 
xcms_peaks_info_less <- xcms_peak_info %>% 
  select(c(Feature, mzmed, mzmin, mzmax, rtmed, rtmin, rtmax))

# Add RAMClustR cluster info 
xcms_peak_info_RC <- xcms_peak_info[RC_res$xcmsOrd,] %>% # fix order 
  mutate(RC = RC_res$featclus)

# Extract xcms peak area/intensity data 
xcms_into <- featureValues(res_cor$PeakID_obj_gr_filled,  value = "into") %>% 
  as.data.frame() %>% 
  rownames_to_column("Feature") # peak area 
xcms_maxint <- featureValues(res_cor$PeakID_obj_gr_filled, value = "maxo") %>% 
  as.data.frame() %>% 
  rownames_to_column("Feature") # max intensity peak

# Prep folders 
dir.create(file.path(out_dir, "Peak_info"))
dir.create(file.path(out_dir, "XCMS_quantification"))
dir.create(file.path(out_dir, "RC_quantification"))

# Save output
write_tsv(xcms_into, file = gzfile(file.path(out_dir, "XCMS_quantification", paste0(file_name,"_XCMS_peak_area.tsv.gz"))))
write_tsv(xcms_maxint, file = gzfile(file.path(out_dir, "XCMS_quantification", paste0(file_name,"_XCMS_peak_max_intensity.tsv.gz"))))
write_tsv(xcms_peak_info, file = gzfile(file.path(out_dir, "Peak_info", paste0(file_name,"_XCMS_peak_info_full.tsv.gz"))))
write_tsv(xcms_peaks_info_less, file = gzfile(file.path(out_dir, "Peak_info", paste0(file_name,"_XCMS_peak_info.tsv.gz"))))
write_tsv(xcms_peak_info_RC, file = gzfile(file.path(out_dir, "Peak_info", paste0(file_name,"_XCMS_peak_info_full_RC.tsv.gz"))))
write_tsv(RC_cluster_res, file = gzfile(file.path(out_dir, "RC_quantification", paste0(file_name,"_RC_peakgr_intensity.tsv.gz"))))
save(xcms_peak_info, file = file.path(out_dir, "Peak_info", "Peak_info_index.RData"))
write_tsv(RC_info, file = file.path(out_dir, "Peak_info", paste0(file_name, "_RC_info.tsv")))

