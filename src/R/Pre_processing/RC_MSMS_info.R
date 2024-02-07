library(tidyverse)
# Load data 
out_dir <- "~/Sara/exAtlas/res/MSMS_info"
MSMS_info <- read.delim("~/Sara/exAtlas/res/MSMS_info/MSMS_summary_info.tsv") 
exAtlas_MS_RC_info <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/Peak_info/exAtlas_MS_RC_info.tsv")
load("~/Sara/exAtlas/res/RAMClustR/MS/exAtlas_MS_RAMClust_res.RData")

# Tag MSMS with 
MSMS_info$ID <- cumsum(!duplicated(MSMS_info[c("precursorScanNum", "precursorMZ")]))
MSMS_info <- MSMS_info %>% 
  dplyr::filter(grepl("QC11", Files)) %>% 
  dplyr::filter(grepl("Iter1", Files)) %>% 
  arrange(seqNum) %>%  
  add_count(ID) %>%
  mutate(newID = ifelse(n == 6, paste(ID, gl(n()/3, 3), sep = "."), ID)) %>% 
  group_by(newID) %>% 
  mutate(precursorScanNum_mz_ID = cur_group_id()) %>% 
  ungroup() %>% 
  dplyr::select(-c(ID, n, newID))
write_tsv(MSMS_info, file.path(out_dir, "MSMS_summary_info_QC11_Iter1_uniqueID.tsv"))

# Make one with all files
MSMS_info2 <- MSMS_info %>% 
  #mutate(retentionTime =  signif(retentionTime, 6),
  #       precursorMZ = signif(precursorMZ, 6)) %>% 
  arrange(desc("seqNum")) %>% 
  group_by(fileIdx, precursorScanNum) %>% 
  mutate(precursorScanNumID = paste(fileIdx, precursorScanNum, gl(n()/3, 3), sep = ".")) %>% 
  as.data.frame()
write_tsv(MSMS_info2, file.path(out_dir, "MSMS_summary_info_uniqueID.tsv"))

# Add topInt 
for(i in 1:nrow(exAtlas_MS_RC_info)){
  my_RC <- exAtlas_MS_RC_info[i,"RC"] %>% as.integer()
  exAtlas_MS_RC_info[i,"mz_topint"] <- RC_res$M.ann[[my_RC]][which.max(RC_res$M.ann[[my_RC]]$int),]$mz
}

# Set parameters 
mz_type = "mz_topint"
ppm_error = 10/1000000
rt_error = 3

# Find matching RC 
MSMS_RC_dat <- lapply(1:nrow(exAtlas_MS_RC_info), function(i){
  # Select mz column to search for productIon 
  my_RC <- exAtlas_MS_RC_info[i,"RC"]
  if(mz_type == "mz"){
    my_mz <- exAtlas_MS_RC_info[i,"mz"] - 1.0073
  } else {
    my_mz <- exAtlas_MS_RC_info[i,"mz_topint"]
  }
  my_rt <- exAtlas_MS_RC_info[i,"rt"]
  
  # Extract Collision Energy levels
  colE <- MSMS_info$collisionEnergy %>% unique()
  
  # Check if we have an MSMS spectra for each energy level 
  res <- lapply(colE, function(b){
    
    # Select a Collision Energy level 
    MSMS_datE <- MSMS_info %>% filter(collisionEnergy == b)
    
    MSMS_datE <- MSMS_datE %>% dplyr::filter(abs(precursorMZ-my_mz)<(ppm_error*my_mz), abs(retentionTime-my_rt)<rt_error)
    
    # Find the closest matching precursorMz 
    my_MSMS <- MSMS_datE[which.min(abs(MSMS_datE$precursorMZ - my_mz)),]
    
    if(!is.null(my_MSMS) & ncol(my_MSMS)>1){
      # Save data 
      my_MSMS %>% mutate(RC = my_RC)
    }
  }) %>% bind_rows()
  res
})
MSMS_RC_dat <- do.call("rbind", purrr::compact(MSMS_RC_dat))

write_tsv(MSMS_RC_dat, file.path(out_dir, "RC_topInt_mz_MSMS_info.tsv"))
#write_tsv(MSMS_RC_dat, file.path(out_dir, "RC_mz_MSMS_info.tsv"))

