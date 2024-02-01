######################################
# exAtlas RC tissue normalization 
######################################
suppressPackageStartupMessages({
  library(tidyverse)
})

# Load data
out_dir <- "~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification"
SPG_int <- read_tsv("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification/SPGs_targeted_RC_peak_intensity.tsv", show_col_types = FALSE)
phenolic_int <- read_tsv("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification/Phenolics/Phenolic_targeted_RC_peak_intensity.tsv", show_col_types = FALSE)
RC_cluster_res <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification/exAtlas_MS_RC_peakgr_intensity.tsv.gz")

file_name <- "exAtlas_MS"
exAtlas_sample_dry_weight <- read_delim("Sara/exAtlas/meta/exAtlas_sample_dry_weight.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
pd <- read_tsv("~/Sara/exAtlas/meta/exAtlas_MS_meta_info.tsv", show_col_types = FALSE)

# Add weight to meta data 
pd_mg <- left_join(pd, exAtlas_sample_dry_weight %>% select(c(Tube_nr, Weight_mg)), by = c(Tube.nr = "Tube_nr")) %>% 
  select(c(Sample, Weight_mg))

# Prep data frame 
SPG_int2 <- SPG_int %>% 
  mutate("Trait_RC" = paste0(Trait, "_RC", RC)) %>% # need unique rownames 
  column_to_rownames("Trait_RC") %>%  
  select(-c(Trait, RC, mz, rt, pre_ion, pre_ion_type)) %>% 
  t() %>% 
  as.data.frame() %>% 
  select_if(~sum(!is.na(.)) > 0) %>% # Remove rows that are completely empty
  rownames_to_column("Sample")

phenol_int2 <- phenolic_int %>% 
  mutate("Trait_RC" = paste0(Trait, "_RC", RC)) %>% # need unique rownames 
  column_to_rownames("Trait_RC") %>%  
  select(-c(Trait, RC, mz, rt, pre_ion, pre_ion_type)) %>% 
  t() %>% 
  as.data.frame() %>% 
  select_if(~sum(!is.na(.)) > 0) %>% # Remove rows that are completely empty
  rownames_to_column("Sample")

untarg_int2 <- RC_cluster_res %>% 
  mutate("Trait_RC" = paste0("RC", RC, "_", mz, "_", rt)) %>% # need unique rownames 
  column_to_rownames("Trait_RC") %>%  
  select(-c(RC, mz, rt, pre_ion, pre_ion_type)) %>% 
  t() %>% 
  as.data.frame() %>% 
  select_if(~sum(!is.na(.)) > 0) %>% # Remove rows that are completely empty
  rownames_to_column("Sample")

# Add weight column 
SPG_int3 <- left_join(pd_mg, SPG_int2, by = "Sample") %>% 
  column_to_rownames("Sample")

phenol_int3 <- left_join(pd_mg, phenol_int2, by = "Sample") %>% 
  column_to_rownames("Sample")

untarg_int3 <- left_join(pd_mg, untarg_int2, by = "Sample") %>% 
  column_to_rownames("Sample")

# Normalize based on weight 
SPG_int_mg <- mutate(SPG_int3,
                     across(.cols = c(!Weight_mg),
                            .fns = ~ . / Weight_mg)) %>% 
  filter_all(any_vars(!is.na(.))) # remove samples with only NA (e.g., BLANK)

phenol_int_mg <- mutate(phenol_int3,
                     across(.cols = c(!Weight_mg),
                            .fns = ~ . / Weight_mg)) %>% 
  filter_all(any_vars(!is.na(.))) # remove samples with only NA (e.g., BLANK)

untarg_int_mg <- mutate(untarg_int3,
                     across(.cols = c(!Weight_mg),
                            .fns = ~ . / Weight_mg))
untarg_int_mg2 <- untarg_int_mg[rowSums(is.na(untarg_int_mg)) != ncol(untarg_int_mg), ] # remove samples with only NA (e.g., BLANK)

# Add additional meta columns 
SPG_int_mg2 <- right_join(pd, SPG_int_mg %>% rownames_to_column("Sample"), by = "Sample")
phenol_int_mg2 <- right_join(pd, phenol_int_mg %>% rownames_to_column("Sample"), by = "Sample")
untarg_int_mg3 <- right_join(pd, untarg_int_mg2 %>% rownames_to_column("Sample"), by = "Sample")

# Save data 
write_tsv(SPG_int_mg2, file = file.path(out_dir, paste0(file_name, "_targeted_peak_intensities_tissue_normalized.tsv")))
write_tsv(phenol_int_mg2, file = file.path(out_dir, "Phenolics", "Phenolic_targeted_peak_intensities_tissue_normalized.tsv"))
write_tsv(untarg_int_mg3, file = file.path(out_dir, paste0(file_name, "_untargeted_peak_intensities_tissue_normalized.tsv")))








