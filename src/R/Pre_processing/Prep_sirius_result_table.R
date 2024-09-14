suppressPackageStartupMessages({
  library(tidyverse)
})

out_dir <- "~/Sara/exAtlas/res/Preprocessing_res/MS/Peak_info/sirius_RC"
compound_identifications <- read.delim("~/Git/exAtlas/src/doc/my_sum5/compound_identifications.tsv")
formula_identifications <- read.delim("~/Git/exAtlas/src/doc/my_sum5/formula_identifications.tsv") %>% 
  mutate(ionMass = round(ionMass,5), 
         retentionTimeInSeconds = round(retentionTimeInSeconds,5))
canopus_formula_summary <- read.delim("~/Git/exAtlas/src/doc/my_sum5/canopus_formula_summary.tsv")

# MSMS data 
MSMS_info <- read.delim("~/Sara/exAtlas/res/MSMS_info/MSMS_summary_info_uniqueID.tsv") %>% 
  dplyr::filter(grepl("QC11", Files)) %>% 
  dplyr::filter(grepl("Iter1", Files)) %>% 
  mutate(precursorMZ = round(precursorMZ,5), 
         retentionTime = round(retentionTime,5))

# Combine sirius formulas with MSMS data
formula_identification_MSMS <- left_join(formula_identifications, MSMS_info, by = c(ionMass = "precursorMZ", retentionTimeInSeconds = "retentionTime")) 
stopifnot(length(unique(formula_identification_MSMS$id)) == length(unique(formula_identifications$id)))

# Add compound identification info to formula MSMS
form_compound_identifications_MSMS <- left_join(formula_identification_MSMS, compound_identifications %>% 
                                                  dplyr::select(-c(rank, SiriusScore, molecularFormula,adduct, ionMass, retentionTimeInSeconds)), by = "id")

# Add compound identification info to formula MSMS
form_cano_compound_identifications_MSMS <- left_join(form_compound_identifications_MSMS, canopus_formula_summary %>% 
                                                  dplyr::select(-c(molecularFormula,adduct)), by = "id")

# Save table 
write_tsv(form_cano_compound_identifications_MSMS, file.path(out_dir, "Sirius_candidates_summary.tsv"))


# Convert Sirius result to RC 
Sirius_candidates_summary <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/Peak_info/sirius_RC/Sirius_candidates_summary.tsv")
load("~/Sara/exAtlas/res/RAMClustR/MS/exAtlas_MS_RAMClust_res.RData")
xcms_peak_info_RC <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/Peak_info/exAtlas_MS_XCMS_peak_info_full_RC_singletons_rm.tsv")

# Set error 
ppm_error <- 10/1000000
rt_error <- 6

# Extract MSMS RC 
XCMS_RC_res <- lapply(1:nrow(Sirius_candidates_summary), function(i){
  # Trait info 
  my_mass <- Sirius_candidates_summary[i,"ionMass"]
  my_rt <- Sirius_candidates_summary[i,"retentionTimeInSeconds"]
  
  # Extract results
  my_res <- xcms_peak_info_RC %>% 
    dplyr::filter(abs(mzmed-my_mass)<(my_mass*ppm_error), abs(rtmed-my_rt)<rt_error) %>% 
    as.data.frame() 
  
  # Combine results 
  xcms_RC_sum <- my_res %>% select(c(Feature, mzmed, mzmin, mzmax, rtmed, rtmin, rtmax, RC))
  
  # Add RC cluster annotations 
  if(nrow(xcms_RC_sum) > 0){
    for(k in 1:nrow(xcms_RC_sum)){
      my_RC <- xcms_RC_sum[k,"RC"]
      xcms_RC_sum[k,"Molecular_ion"] <- RC_res$M[my_RC]
      xcms_RC_sum[k,"Precursor_ion"] <- RC_res$precursor.mz[my_RC]
      xcms_RC_sum[k,"Precursor_type"] <- RC_res$precursor.type[my_RC]
      xcms_RC_sum[k,"RC_rt"] <- RC_res$clrt[my_RC]
      xcms_RC_sum[k,"Trait"] <-  Sirius_candidates_summary[i,"name"]
      xcms_RC_sum[k,"Targ_mz"] <- my_mass
      xcms_RC_sum[k,"Targ_rt"] <- my_rt
      xcms_RC_sum[k,"precursorScanNumID"] <- Sirius_candidates_summary[i,"precursorScanNumID"]      
      xcms_RC_sum[k,"id"] <- Sirius_candidates_summary[i,"id"]
      xcms_RC_sum[k,"collisionE"] <- Sirius_candidates_summary[i,"collisionEnergy"] 
      xcms_RC_sum[k,"smiles"] <- Sirius_candidates_summary[i,"smiles"]
      xcms_RC_sum[k,"NPC.superclass"] <- Sirius_candidates_summary[i,"NPC.superclass"]
      xcms_RC_sum[k,"NPC.superclass.Probability"] <- Sirius_candidates_summary[i,"NPC.superclass.Probability"]
    }
  }
  xcms_RC_sum
}) %>% bind_rows()

write_tsv(XCMS_RC_res, file.path(out_dir, "Sirius_candidates2RC_singletons_rm.tsv"))


