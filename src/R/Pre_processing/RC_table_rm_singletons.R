suppressPackageStartupMessages({
  library(tidyverse)
})

# Load data 
untarg_int <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification/exAtlas_MS_untargeted_peak_intensities_tissue_normalized.tsv")
peak_info_RC <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/Peak_info/exAtlas_MS_XCMS_peak_info_full_RC.tsv.gz")
load("~/Sara/exAtlas/res/RAMClustR/MS/exAtlas_MS_RAMClust_res.RData")
RC_cluster_res <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification/exAtlas_MS_RC_peakgr_intensity.tsv.gz")

# Identify columns starting with "RC"
rc_columns <- grep("^RC", names(untarg_int), value = TRUE)

# Modify column names
new_names <- sub("_.*", "", rc_columns)

# Rename the columns in the dataframe
names(untarg_int)[names(untarg_int) %in% rc_columns] <- new_names

# Find singleton boarder
stopifnot(nrow(RC_res$M.ann[[4125]]) > 1)
stopifnot(nrow(RC_res$M.ann[[4126]]) > 1)

# Remove singletons 
untarg_int_sel <- untarg_int %>% select(c(Sample:Weight_mg, paste0("RC", 1:4125))) 

# Cut RT 
RC_keep <- peak_info_RC %>% filter(RC %in% str_remove(grep("^RC", names(untarg_int_sel), value = TRUE), "RC")) %>% 
  filter(rtmax < 540, rtmin > 30) %>% pull(RC) %>% unique() %>% sort()

# Save final table 
untarg_int_final <- untarg_int_sel %>% select(c(Sample:Weight_mg, paste0("RC", RC_keep))) 
write_tsv(untarg_int_final, file.path("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification", "exAtlas_MS_untargeted_peak_intensities_tissue_normalized_singletons_rm.tsv"))

peak_info_RC_final <- peak_info_RC %>% filter(RC %in% RC_keep)
write_tsv(peak_info_RC_final, file.path("~/Sara/exAtlas/res/Preprocessing_res/MS/Peak_info", "exAtlas_MS_XCMS_peak_info_full_RC_singletons_rm.tsv"))

RC_cluster_res_final <- RC_cluster_res %>% filter(RC %in% RC_keep)
write_tsv(RC_cluster_res_final, file.path("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification", "exAtlas_MS_RC_peakgr_intensity_singletons_rm.tsv"))
