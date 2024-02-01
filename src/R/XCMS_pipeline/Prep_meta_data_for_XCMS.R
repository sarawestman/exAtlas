suppressPackageStartupMessages({
  require(readr)
  require(tidyverse)
  require(dplyr)
})

## Prep pd dataframe for XCMS ##
out_dir <- "~/Sara/exAtlas/meta"
file_name <- "exAtlas_MS_meta_info"
meta_info <- read.csv("~/Sara/exAtlas/meta/exAtlas_tube_labels_R_fixed.csv", sep=";")
Worklist <- read.csv("~/Sara/exAtlas/meta/Worklist.csv", sep=";")
sub.folders <- list.dirs("~/Sara/exAtlas/data/MS", recursive = TRUE)[-1]
fls <- list.files(path = sub.folders, pattern = "*.mzML.gz", full.names = TRUE, recursive = FALSE)

# sanity check 
stopifnot(str_remove(basename(fls), ".mzML.gz") %in% Worklist$Sample)

# Create meta file 
pd_inc <- Worklist[Worklist$Sample %in% str_remove(basename(fls), ".mzML.gz"),] 
pd_combo <- left_join(pd_inc %>% mutate(Sample_caps = toupper(Sample)), meta_info %>% 
                        mutate(Sample_caps = toupper(Name_NEG)), by = "Sample_caps") %>% select(-c(Name_NEG,Mode))
pd_fix <- pd_combo %>% mutate(Tissue = ifelse(is.na(Tissue), Type, Tissue),
                              Tissue_exp = ifelse(is.na(Tissue_exp), Type, Tissue_exp),
                              Experiment = ifelse(is.na(Experiment), Type, Experiment))
pd_fix_ord <- pd_fix[order(match(pd_fix$Sample,str_remove(basename(fls), ".mzML.gz"))),]
stopifnot(pd_fix_ord$Sample == str_remove(basename(fls), ".mzML.gz"))
pd <- pd_fix_ord %>% mutate(File = fls, Genotype_nr = paste0(Tissue_exp, "_",Genotype)) %>% arrange(Injection) %>% dplyr::select(Sample, everything())


# Check - For upcoming scripts 
if(!"File" %in% colnames(pd) | !"Genotype" %in% colnames(pd) | !"Genotype_nr" %in% colnames(pd) | !"Injection" %in% colnames(pd) | !"Sample" %in% colnames(pd)){
  stop("Column(s) missing in pd, fix column(s) or adapt script")
}

# Save data 
write_tsv(pd, file = file.path(out_dir, paste0(file_name,".tsv")))
