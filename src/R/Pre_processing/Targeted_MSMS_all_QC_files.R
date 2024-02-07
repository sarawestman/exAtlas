suppressPackageStartupMessages({
  library(tidyverse)
  library(MSnbase)
  library(ggtext)
  source("~/Git/exAtlas/src/R/CSPP/Explore2/MSMS_functions.R")
})

# Set parameters 
# 0.2 or even 0.5 ppm
ppm_error = 10/1000000
rt_error = 5
nlabel = 40 
out_dir <- "~/Sara/exAtlas/res/Preprocessing_res/MS/Plots/MSMS_targeted/allFiles/SPGs/RT_5_ppm_10/smaller"

targ_cpd <- read.delim("~/Sara/exAtlas/meta/Targeted_metabolites/SPG_targeted_complete_list_adducts.tsv")
MSMS_info <- read.delim("~/Sara/exAtlas/res/MSMS_info/MSMS_summary_info_uniqueID.tsv") 
  
# Find targeted FA 
MSMS_dat_targ_FA <- lapply(1:nrow(targ_cpd), function(i){
  my_mz <- targ_cpd[i,"Mass_FA"] 
  my_rt <- targ_cpd[i,"RT"]
  
  form_id_sel <- MSMS_info %>% dplyr::filter(abs(precursorMZ-my_mz)<(ppm_error*my_mz), abs(retentionTime-my_rt)<rt_error)

  # Find the closest matching precursorMz 
  form_id_sel2 <- form_id_sel %>% dplyr::arrange(desc(precursorIntensity)) %>% slice(1)
  
  if(nrow(form_id_sel2)>0){
    form_id_sel2 %>% mutate(Trait = targ_cpd[i,"Cpd"], 
                            Type = "FA", 
                            Formula = targ_cpd[i,"Formula"], 
                            Origin = targ_cpd[i,"Comments"])
  }
}) %>% bind_rows() %>% distinct()

# Find targeted H
MSMS_dat_targ_H <- lapply(1:nrow(targ_cpd), function(i){
  my_mz <- targ_cpd[i,"Mass_H"]
  my_rt <- targ_cpd[i,"RT"]
  
  form_id_sel <- MSMS_info %>% dplyr::filter(abs(precursorMZ-my_mz)<(ppm_error*my_mz), abs(retentionTime-my_rt)<rt_error)

  # Find the closest matching precursorMz 
  form_id_sel2 <- form_id_sel %>% dplyr::arrange(desc(precursorIntensity)) %>% slice(1)
  
  if(nrow(form_id_sel2)>0){
    form_id_sel2 %>% mutate(Trait = targ_cpd[i,"Cpd"], 
                            Type = "H", 
                            Formula = targ_cpd[i,"Formula"], 
                            Origin = targ_cpd[i,"Comments"]) 
  }
}) %>% bind_rows() %>% distinct()

# Combine tables
(MSMS_targ <- rbind(MSMS_dat_targ_FA, MSMS_dat_targ_H))
write_tsv(MSMS_targ, file.path(out_dir, "Targeted_MSMS_info.tsv"))

# Plot targeted MSMS - including all collision energies 
for(k in 1:nrow(MSMS_targ)){
  my_id <- MSMS_targ[k,"precursorScanNumID"]
  MSMS_sel <- MSMS_info %>% dplyr::filter(precursorScanNumID == my_id) %>% mutate(Trait = MSMS_targ[k,"Trait"], Type =  MSMS_targ[k,"Type"])
  plot_list <- list()
  
  # If we have MSMS spectra for a Substrate, plot the spectra
  if(nrow(MSMS_sel) > 0){
    
    # Store plot in a list 
    for(i in 1:nrow(MSMS_sel)){
      
      # Extract SeqNum ID 
      my_seqN <- MSMS_sel[i, "seqNum"]
      
      # Read SeqNum MSMS data 
      my_MSMS <- readMSData(MSMS_sel[i,"Files"], mode = "onDisk")
      
      # Plot the spectra
      my_mz <- my_MSMS[[my_seqN]]@mz
      txt_lab <- round(my_MSMS[[my_seqN]]@mz,2)
      
      y_int <- my_MSMS[[my_seqN]]@intensity
      
      x_lab = my_MSMS[[my_seqN]]@mz[order(my_MSMS[[my_seqN]]@intensity, decreasing = TRUE)[1:nlabel]]
      y_lab = my_MSMS[[my_seqN]]@intensity[order(my_MSMS[[my_seqN]]@intensity, decreasing = TRUE)[1:nlabel]]
      txt_lab2 = round(my_MSMS[[my_seqN]]@mz[order(my_MSMS[[my_seqN]]@intensity, decreasing = TRUE)[1:nlabel]],2)
      
      # Plot 
      my_plot <- MSnbase::plot(my_MSMS[[my_seqN]]) + 
        annotate('richtext', 
                 x = x_lab, 
                 y = y_lab, 
                 label = txt_lab2, 
                 fill = NA, 
                 label.color = NA) + 
        theme_classic() + 
        labs(title = paste0("Name: ", paste0(MSMS_sel[i, "Trait"], " (",paste0(MSMS_sel[i, "Type"]), ") "), "File: ", basename(MSMS_sel[i, "Files"])),
             subtitle = paste0("MZ: ", round(MSMS_sel[i, "precursorMZ"], 2), ". RT: ", round(MSMS_sel[i, "retentionTime"], 2),". Collision E: ", MSMS_sel[i,"collisionEnergy"]))
      plot_list[[paste0(my_id, "_", i)]] <- print(my_plot)
    }
    
    message("Plotting") # 15, 30 vs 30, 50 
    pdf(file.path(out_dir, paste0(my_id,".pdf")), width = 50, height = 30)
    for(i in sub("_[^_]+$", "", names(plot_list)) %>% unique()){
      cnames <- names(plot_list)[str_detect(names(plot_list), i)]
      nplots <- names(plot_list)[str_detect(names(plot_list), i)] %>% length()
      
      # Adapt plot columns depending on number of Collision Energy levels 
      if(nplots == 1){
        gridExtra::grid.arrange(plot_list[cnames][[1]], nrow=3)
      }
      if(nplots == 2){
        gridExtra::grid.arrange(plot_list[cnames][[1]], plot_list[cnames][[2]], nrow=3) 
      }
      if(nplots == 3){
        gridExtra::grid.arrange(plot_list[cnames][[1]], plot_list[cnames][[2]], plot_list[cnames][[3]], nrow=3)   
      }
    }
    graphics.off()
  }
}

# Sanity check 
for(k in 1:nrow(MSMS_targ)){
  my_id <- MSMS_targ[k,"precursorScanNumID"]
  MSMS_sel <- MSMS_info %>% dplyr::filter(precursorScanNumID == my_id) %>% mutate(Trait = MSMS_targ[k,"Trait"], Type =  MSMS_targ[k,"Type"])
  message(paste("Nr of precursorScanNumID found:", nrow(MSMS_sel)))
}
