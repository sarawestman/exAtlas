suppressPackageStartupMessages({
  library(tidyverse)
  library(MSnbase)
})

out_dir <- "~/Sara/exAtlas/res/MSMS_info"

# Read files
MSMS_files <- data.frame(Files = list.files(c("~/Sara/exAtlas/data/MSMS/B5", 
                                               "~/Sara/exAtlas/data/MSMS/B1", 
                                               "~/Sara/exAtlas/data/MSMS/mzML_Tissue_MSMS"
                                               ), full.names = TRUE)) %>% rowid_to_column() %>% filter(!grepl("Blank", Files)) %>% mutate(rowid = seq(1:nrow(.)))

# Read MSMS data
data_MSMS <- readMSData(MSMS_files$Files, mode = "onDisk")
MSMS_precursor <- left_join(fData(data_MSMS), MSMS_files, by = c(fileIdx = "rowid")) %>% 
  select(fileIdx, 
         seqNum, 
         msLevel, 
         precursorScanNum, 
         collisionEnergy, 
         precursorMZ, 
         precursorIntensity, 
         retentionTime, 
         Files) 
MSMS_info <- MSMS_precursor[!is.na(MSMS_precursor$precursorMZ),]

# Save data 
write_tsv(MSMS_info, file.path(out_dir, "MSMS_summary_info.tsv"))