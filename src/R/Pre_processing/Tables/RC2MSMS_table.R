suppressPackageStartupMessages({
  library(tidyverse)
  library(MSnbase)
  library(ggtext)
  source("~/Git/exAtlas/src/R/CSPP/Explore2/MSMS_functions.R")
})

# Load data 
out_dir <- "~/Sara/exAtlas/res/MSMS_info"
MSMS_info <- read.delim("~/Sara/exAtlas/res/MSMS_info/MSMS_summary_info_uniqueID.tsv")
xcms_peak_info_RC <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/Peak_info/exAtlas_MS_XCMS_peak_info_full_RC_singletons_rm.tsv")
load("~/Sara/exAtlas/res/RAMClustR/MS/exAtlas_MS_RAMClust_res.RData")
sirius_RC_MSMS <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/Peak_info/sirius_RC/Sirius_candidates2RC_singletons_rm.tsv") %>%
  filter(NPC.superclass.Probability > 0.8) %>% 
  group_by(RC) %>% 
  top_n(1, NPC.superclass.Probability) %>% # The NPC superclass with highest probability score will be used to represent the whole RC cluster. 
  ungroup()

# Add additional RC info to table 
xcms_peak_info_RC2 <- lapply(unique(xcms_peak_info_RC$RC), function(x){
  RC_info <- RC_res$M.ann[[x]]
  names(RC_info) <- paste0("RC_", names(RC_info))
  xcms_info <- xcms_peak_info_RC %>% filter(RC == x) 
  
  if(nrow(xcms_info) > 1){
  xcms_info <- lapply(1:nrow(xcms_info), function(i){
    my_mz <- xcms_info[i, "mzmed"]
    frag_index <- which.min(abs(RC_info$RC_mz - my_mz))
    RC_add <- RC_info[frag_index,]
    xcms_info[i,] %>% bind_cols(RC_add)
    }) %>% bind_rows()
  }
  xcms_info
  }) %>% bind_rows()

# Set errors 
ppm_error = 10/1000000
rt_error = 5

# Find MSMS spectra of xcms feature
MSMS_xcms_dat <- lapply(1:nrow(xcms_peak_info_RC2), function(i){
  my_mz <- xcms_peak_info_RC2[i,"mzmed"] 
  my_rt <- xcms_peak_info_RC2[i,"rtmed"]

  form_id_sel <- MSMS_info %>% dplyr::filter(abs(precursorMZ-my_mz)<(ppm_error*my_mz), abs(retentionTime-my_rt)<rt_error)
  
  if(nrow(form_id_sel) > 0){
    form_id_sel$rt_tmp <- abs(form_id_sel$retentionTime-my_rt)
    form_id_sel$mz_tmp <- abs(form_id_sel$precursorMZ-my_mz)
  
    # Find matching precursorMz - select the one with highest intensity (that still passed the set thresholds)
    form_id_sel2 <- form_id_sel %>% dplyr::arrange(desc(precursorIntensity)) %>% head(n = 1)

    if(nrow(form_id_sel2) > 0){
      form_id_sel2 %>% mutate(Feature = xcms_peak_info_RC2[i,"Feature"]) # Add feature ID 
    }
    }
  }) %>% 
  bind_rows() %>% 
  distinct() %>% 
  mutate(dup = duplicated(precursorScanNumID)) %>% 
  group_by(precursorScanNumID) %>% 
  slice(which.min(rt_tmp)) %>% # If MSMS spectra is matched with multiple features, select the one with the closest matching RT   
  select(-c(dup))

# Add MSMS info to xcms peak table 
MSMS_xcms_final <- left_join(xcms_peak_info_RC2, MSMS_xcms_dat, by = "Feature") %>% select(c(RC:RC_label, ms_level, Feature:peakidx, msLevel, fileIdx, seqNum, precursorScanNum:mz_tmp))

# Add info from other collision levels
colnames(MSMS_info) <- paste0("MSMS_", colnames(MSMS_info))
MSMS_xcms_final2 <- left_join(MSMS_xcms_final %>% select(-c(peakidx, collisionEnergy, msLevel, seqNum, retentionTime, RC_label, Files, rt_tmp, mz_tmp, fileIdx:precursorIntensity)), MSMS_info, by = c("precursorScanNumID" = "MSMS_precursorScanNumID")) 

# Extract Features with MSMS info
MSMS_hits <- MSMS_xcms_final2[!is.na(MSMS_xcms_final2$MSMS_Files),] 
perc_excl = 90 # Set an intensity threshold in percentage (will remove X % of the lowest intensity ions for an individual feature)

# Add MSMS product ion m/z and int data to features
MSMS_xcms_final3 <- lapply(1:nrow(MSMS_hits), function(i){
  
  # Extract SeqNum ID 
  my_seqN <- MSMS_hits[i, "MSMS_seqNum"]

  # Read SeqNum MSMS data 
  my_MSMS <- readMSData(MSMS_hits[i,"MSMS_Files"], mode = "onDisk")

  # Plot the spectra
  my_mz <- my_MSMS[[my_seqN]]@mz
  my_int <- my_MSMS[[my_seqN]]@intensity
  
  # Store min and max int info 
  maxInt <- max(my_int)
  minInt <- min(my_int)
  
  # Filter mz and int based on intensity threshold
  nInc <- round(((100-perc_excl)/100)*length(my_int), 0)
  int_lab = round(my_int[order(my_int, decreasing = TRUE)[1:nInc]], 4)
  mz_lab = round(my_mz[order(my_int, decreasing = TRUE)[1:nInc]], 4)
  IntThr <- min(int_lab)
  
  MSMS_hits[i,] %>% select(c(Feature, precursorScanNumID, MSMS_collisionEnergy)) %>% 
    mutate(MSMS_product_mz = paste(mz_lab,  collapse = ';'), 
           MSMS_product_int = paste(int_lab,  collapse = ';'),
           MSMS_int_threshold = IntThr,
           MSMS_minInt = minInt,
           MSMS_maxInt = maxInt,
           MSMS_nr_frag_inc = nInc,
           MSMS_nr_frag_rm = length(my_int)-length(int_lab)
           )
  }) %>% bind_rows()

MSMS_xcms_final4 <- left_join(MSMS_xcms_final2, MSMS_xcms_final3, by = c("Feature", "precursorScanNumID", "MSMS_collisionEnergy"))

# Add sirius annotation
MSMS_xcms_sirius <- left_join(MSMS_xcms_final4, sirius_RC_MSMS %>% dplyr::rename("sirius_precursorScanNumID" = precursorScanNumID) %>% 
                                select(c(Feature, id, sirius_precursorScanNumID, collisionE, smiles:NPC.superclass.Probability)) %>% 
                                dplyr::rename("sirius_id" = "id", 
                                              "sirius_smile" = "smiles", 
                                              "sirius_NPC_superclass" = "NPC.superclass", 
                                              "sirius_NPC_superclass_probability" = "NPC.superclass.Probability"), 
                              by = c("Feature" = "Feature", 
                                     "MSMS_collisionEnergy" = "collisionE")) %>% 
  dplyr::rename("RC_id" = "RC", "MSMS_path" = "MSMS_Files") %>% 
  mutate(RC_id = paste0("RC", RC_id)) %>% 
  select(c(RC_id:Stem_In_vitro_Elite, MSMS_msLevel, precursorScanNumID:MSMS_seqNum, MSMS_precursorScanNum:sirius_NPC_superclass_probability))

# Save results
write_tsv(MSMS_xcms_sirius, file.path(out_dir, paste0("exAtlas_MS_XCMS_peak_info_full_RC_MSMS_singletons_rm_sirius_MSMSint_", perc_excl, "_percent_excluded.tsv.gz")))

# Later column rearrangements 
MSMS_xcms_sirius <- read.delim("~/Sara/exAtlas/res/MSMS_info/exAtlas_MS_XCMS_peak_info_full_RC_MSMS_singletons_rm_sirius_MSMSint_90_percent_excluded.tsv.gz") %>% mutate(MSMS_file = basename(MSMS_path))
write_tsv(MSMS_xcms_sirius %>% select(c(RC_id:MSMS_path, MSMS_file, MSMS_product_mz:sirius_NPC_superclass_probability)), file = "~/Sara/exAtlas/res/MSMS_info/exAtlas_MS_XCMS_peak_info_full_RC_MSMS_singletons_rm_sirius_MSMSint_90_percent_excluded.tsv.gz")


