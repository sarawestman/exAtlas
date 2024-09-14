suppressPackageStartupMessages({
  library(tidyverse)
})

# Load data 
out_dir <- "~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification/Proanthocyanidins"
xcms_peak_info_RC <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/Peak_info/exAtlas_MS_XCMS_peak_info_full_RC_singletons_rm.tsv")
load("~/Sara/exAtlas/res/RAMClustR/MS/exAtlas_MS_RAMClust_res.RData")
RC_cluster_res <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification/exAtlas_MS_RC_peakgr_intensity_singletons_rm.tsv")

# Include only SPGs with validated spectra 
#targ_df <- read.delim("~/Sara/exAtlas/meta/Targeted_metabolites/SPG_targeted_complete_list_adducts.tsv")[c(1,3,7,8,4,14,26,19,31,2,22,21,6,11,12,16,17,24),]
#targ_df <- read.delim("~/Sara/exAtlas/meta/Targeted_metabolites/Phenolic_compounds_adducts.tsv")
targ_df <- read.delim("~/Sara/exAtlas/meta/Targeted_metabolites/Proanthocyanidins_targeted_adducts.tsv")
rules <- read.delim("~/Sara/exAtlas/meta/Adduct_ion_rules.tsv")

# Set parameters 
targeted_name <- "Proanthocyanidins"
ppm_error <- 10/1000000
rt_error <- 6

# Save targeted param 
targ_param <- list("ppm_error" = ppm_error,
                   "RT_error_sec" = rt_error)

# Extract targeted data 
XCMS_RC_res1 <- lapply(1:nrow(targ_df), function(i){
  # Trait info 
  my_mass_H <- targ_df[i,"Mass_H"]
  my_mass_FA <- targ_df[i,"Mass_FA"]
  my_rt <- targ_df[i,"RT"]
  my_trait <- targ_df[i,"Cpd"]
  
  # Extract results
  # H
  my_res_H <- xcms_peak_info_RC %>% 
    mutate(Type = "H") %>% 
    dplyr::filter(abs(mzmed-my_mass_H)<(my_mass_H*ppm_error), abs(rtmed-my_rt)<rt_error) %>% 
    as.data.frame() 
  
  # FA
  my_res_FA <- xcms_peak_info_RC %>% 
    mutate(Type = "FA") %>% 
    dplyr::filter(abs(mzmed-my_mass_FA)<(my_mass_FA*ppm_error), abs(rtmed-my_rt)<rt_error) %>%
    as.data.frame() 
  
  # Combine results 
  xcms_RC_sum <- rbind(my_res_H, my_res_FA) %>% 
    select(c(Feature, mzmed, mzmin, mzmax, rtmed, rtmin, rtmax, RC, Type))
  
  # Add RC cluster annotations 
  if(nrow(xcms_RC_sum) > 0){
    for(k in 1:nrow(xcms_RC_sum)){
      my_RC <- xcms_RC_sum[k,"RC"]
      xcms_RC_sum[k,"Molecular_ion"] <- RC_res$M[my_RC]
      xcms_RC_sum[k,"Precursor_ion"] <- RC_res$precursor.mz[my_RC]
      xcms_RC_sum[k,"Precursor_type"] <- RC_res$precursor.type[my_RC]
      xcms_RC_sum[k,"Trait"] <- my_trait
    }
  }
  if(nrow(xcms_RC_sum) > 0){
  xcms_RC_sum
  }
}) 

# Quality check data 
qc_traits <- lapply(XCMS_RC_res1, function(x){
  if(!is.null(x)){
    my_RC <- x %>% pull(RC) %>% unique()
    my_trait <- x %>% pull(Trait)
    
    if(length(my_RC) > 1){
      my_trait <- x %>% pull(Trait) %>% unique()
      message(paste(my_trait, "has multiple RC's"))
      x %>% mutate(Note = "MultiRC")
      } else if(nrow(x) == 1){
        my_trait <- x %>% pull(Trait)
        message(paste(my_trait, "is missing one adduct type"))
        x %>% mutate(Note = "SingleRC")
    }
  }
})

XCMS_RC_res <- do.call(rbind, XCMS_RC_res1) %>% 
  dplyr::select(c(Trait, RC, Feature, everything())) 

XCMS_RC_res_targ <- left_join(XCMS_RC_res, targ_df, by = c(Trait = "Cpd")) 

XCMS_RC_res_targ_diff <- XCMS_RC_res_targ %>%
  mutate(Mass_diff = case_when(Type == "FA" ~ abs(Mass_FA - mzmed),
                               Type == "H" ~ abs(Mass_H - mzmed)),
         RT_diff = abs(RT - rtmed))

#########################################
# Targeted: Peak cluster intensity data 
#########################################
RC_trait <- XCMS_RC_res %>% 
  select(c(Trait, RC)) %>% 
  distinct()
targeted_data <- left_join(RC_trait, RC_cluster_res, by = "RC")

##############################
# Explore RC cluster content  
#############################
# Extract targeted RC clusters - check which peaks they contain 
RC_targ_peaks <- lapply(1:nrow(RC_trait), function(i){
  my_RC <- RC_trait[i,"RC"]
  my_trait <- RC_trait[i,"Trait"]
  
  if(my_RC != 0){
      my_res <- xcms_peak_info_RC %>% filter(RC == my_RC) %>% 
        dplyr::select(c(Feature, mzmed, mzmin, mzmax, rtmed, rtmin, rtmax, RC)) %>% 
        mutate(Trait = my_trait) %>% 
        arrange(mzmed) 
      my_res
      }
})

# Quick version: Add putative annotations (this one is likely to be incomplete - only one column - annotations might overwrite each other)
RC_targ_peaks_anno_short <- lapply(RC_targ_peaks, function(k){
  if(!is.null(k)){
    for(i in 1:nrow(rules)){
      # Extract all mzmed in a cluster 
      mzmed_gr <- k$mzmed
      
      # Extract a rule 
      mzdif <- rules[i,"massdiff"]
      my_mult <- rules[i,"Mult"]
      
      # Calculate mass of rule compound
      for(j in 1:nrow(k)){
        mzmed_rule <- k[j,"mzmed"]*my_mult+mzdif # IMPORTANT! nmol and charge effects on mz 
      
        # Check if rule compound have a match with a mzmed cluster compound
        for(t in mzmed_gr){
          if(abs(mzmed_rule-t)<(t*ppm_error)){
            k[j,"M_mzdif"] <- mzdif
            k[j,"M_derived"] <- t
            k[j,"M_rule"] <-  rules[i,"name"]
          }
        }
        }
      }
  k
  }
  })

# Full version: Add putative annotations (this one is complete - no overwriting each other)
RC_targ_peaks_anno_complete <- lapply(RC_targ_peaks, function(k){
    if(!is.null(k)){
      for(i in 1:nrow(rules)){
        # Extract all mzmed in a cluster 
        mzmed_gr <- k$mzmed
        
        # Extract a rule 
        mzdif <- rules[i,"massdiff"]
        my_mult <- rules[i,"Mult"]
        
        # Calculate mass of rule compound
        for(j in 1:nrow(k)){
          mzmed_rule <- k[j,"mzmed"]*my_mult+mzdif
          
          # Check if rule compound have a match with a mzmed cluster compound
          for(t in mzmed_gr){
            if(abs(mzmed_rule-t)<(t*ppm_error)){
              k[j,paste0(rules[i,"name"])] <- t 
            }
          }
        }
      }
      k 
    }
})

# Save data 
save(targ_param, XCMS_RC_res, targ_df, targeted_data, XCMS_RC_res_targ_diff, RC_targ_peaks, RC_targ_peaks_anno_short, RC_targ_peaks_anno_complete, file = file.path(out_dir, paste0(targeted_name, "_targeted_RC_putative_annotations.RData")))
write_tsv(targeted_data, file = file.path(out_dir, paste0(targeted_name, "_targeted_RC_peak_intensity.tsv")))

