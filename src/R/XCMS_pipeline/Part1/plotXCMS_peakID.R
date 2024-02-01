#########################################################
# XCMS - Part 1: Plot XCMS peak identification results
#########################################################
source("~/Git/exAtlas/src/R/XCMS_pipeline/XCMS_functions.R")
suppressPackageStartupMessages({
  library(xcms)
  library(magrittr)
})
tictoc::tic()

################################################           
# Load files  
################################################
# Output directory
out_dir <- "~/Sara/exAtlas/res/XCMS/MS/Peak_identification"
in_file <- "~/Sara/exAtlas/res/XCMS/MS/Peak_identification/exAtlas_MS_XCMS_peak_identification.RData"
targ_file <- "~/Sara/exAtlas/meta/Targeted_metabolites/SPG_targeted_complete_list_adducts.tsv"

# Targeted data 
targ_df <- read.delim(targ_file)

# Load xcms data 
load(in_file)

################################################           
# Set XCMS parameters  
################################################
# Which plot do you want to make?
targ_plot = TRUE 
untarg_plot = FALSE
nrandom = 50 # How many untargeted peaks should be plotted?

# Plot results from intensity filtering? 
ReFine = TRUE

# Plot results from alignment? 
Align = FALSE

# Select ppm error for targeted peak identification 
ppm_error = 10/1000000
rt_error = 6

################################################           
# Check assumptions
################################################
# Check that targeted columns exist
if(targ_plot){
  if(!"Mass" %in% colnames(targ_df) | !"Mass_FA" %in% colnames(targ_df) | !"Mass_H" %in% colnames(targ_df) | !"Cpd" %in% colnames(targ_df) | !"RT" %in% colnames(targ_df)){
    stop("Column(s) missing in targ_df, fix column(s) or adapt script")
  }
}

# Check if output directory exists 
if(!dir.exists(out_dir)){
  stop(paste("Output directory does not exist:", out_dir))
} else{
  message(paste("Input file:", in_file))
  message(paste("Targeted file used:", targ_file))
  dir.create(file.path(out_dir, "Plots"), showWarnings = FALSE)
  message(paste("Output directory:", file.path(out_dir, "Plots")))
}

# Set parallelisation - nr of cores
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1) {
  stop("Nr of cores for parallelisation must be supplied", call.=FALSE)
} else {
  nCore <- args[1]
  message("Nr of cores: ", nCore)
}

##################################################
# Plot - Sanity check of peaks (random/targeted)
##################################################
# Peak identification check: Random set of peaks
if(untarg_plot){
  message("Plot untargeted peaks")
  # Pick untargeted peaks 
  peakID_rand <- chromPeaks(res_PeakID$Peak_data) %>% 
    rownames() %>% 
    sample(., nrandom)
  peak_rand_df <- chromPeaks(res_PeakID$Peak_data)[peakID_rand,]
  
  # Plotting 
  pdf(file.path(out_dir, "Plots", "XCMS_untargeted_peaks.pdf"), width = 12, height = 12)
  lapply(1:nrow(peak_rand_df), function(i){
    peak_ID <- rownames(peak_rand_df)[[i]]
    trait <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data)) %>% 
      dplyr::filter(abs(mz-peak_rand_df[i,"mz"])<(peak_rand_df[i,"mz"]*ppm_error), abs(rt-peak_rand_df[i, "rt"])<rt_error) 
    if(nrow(trait) != 0){
      peak_plot(res_PeakID$Peak_data, trait, peak_ID, nCore)
    }
  })
  graphics.off()
}

# Peak identification check: Targeted peaks
if(targ_plot){
  message("Plot targeted peaks")
  dir.create(file.path(out_dir, "Plots", "Targeted_peaks"), showWarnings = FALSE)
  lapply(1:nrow(targ_df), function(i){
    
    # Prepare pdf
    trait_name <- targ_df[i,"Cpd"]
    fname <- file.path(out_dir, "Plots", "Targeted_peaks", paste0("XCMS_", i, "_", trait_name, ".png"))
    png(fname, width = 10, height = 5, units = 'in', res = 600)
    par(mfrow = c(1, 2), mar = c(10, 5, 5, 5))
    
    # Extract ions
    trait_FA <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data)) %>% 
      dplyr::filter(abs(mz-targ_df[i,"Mass_FA"])<(targ_df[i,"Mass_FA"]*ppm_error), abs(rt-targ_df[i, "RT"])<rt_error) 
    trait_H <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data)) %>% 
      dplyr::filter(abs(mz-targ_df[i,"Mass_H"])<(targ_df[i,"Mass_H"]*ppm_error), abs(rt-targ_df[i, "RT"])<rt_error)
    
    # Plot the chromatograms (FA)
    if(nrow(trait_FA) != 0){
      my_title <- paste(trait_name, "[M+HCOO]-")
      peak_plot(res_PeakID$Peak_data, trait_FA, my_title, nCore)
    }
    # Plot the chromatograms (H)
    if(nrow(trait_H) != 0){
      my_title <- paste(trait_name, "[M-H]-")
      peak_plot(res_PeakID$Peak_data, trait_H, my_title, nCore)
    }
    graphics.off()
  })
}

#######################################
# Plot - Sanity check (filtered Peak_data)
#######################################
# Check X random peaks that were removed (refined peak data) 
if(ReFine){
  message("Plot removed peaks")
  # Extract removed peak IDs
  Peak_data_peakID <- chromPeaks(res_PeakID$Peak_data) %>% 
    rownames()
  Peak_data_refined_peakID <- chromPeaks(res_PeakID$Peak_data_refined) %>% 
    rownames()
  inf_rm_peaks <- Peak_data_peakID[!Peak_data_peakID %in% Peak_data_refined_peakID]
  seed <- 12345
  lowint_rand <- sample(inf_rm_peaks, nrandom)
  peak_inf_df <- chromPeaks(res_PeakID$Peak_data)[lowint_rand,]
  
  # Plotting
  pdf(file.path(out_dir, "Plots", "XCMS_removed_peaks.pdf"), width = 12, height = 12)
  lapply(1:nrow(peak_inf_df), function(i){
    peak_ID <- rownames(peak_inf_df)[[i]]
    trait <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data)) %>% 
      dplyr::filter(abs(mz-peak_inf_df[i,"mz"])<(peak_inf_df[i,"mz"]*ppm_error), abs(rt-peak_inf_df[i, "rt"])<rt_error) 
    if(nrow(trait) != 0){
      peak_plot(res_PeakID$Peak_data, trait, paste(peak_ID, "removed"), nCore)
    }
  })
  graphics.off()
  
  # Targeted peaks 
  if(targ_plot){
    message("Plot targeted peaks (refined)")
    dir.create(file.path(out_dir, "Plots", "Refined_targeted_peaks"), showWarnings = FALSE)
    lapply(1:nrow(targ_df), function(i){
      
      # Prepare pdf 
      trait_name <- targ_df[i,"Cpd"]
      fname <- file.path(out_dir, "Plots", "Refined_targeted_peaks", paste0("Refined_", i, "_", trait_name, ".png"))
      png(fname, width = 10, height = 4, units = 'in', res = 600)
      par(mfrow = c(1, 2))
      
      # Extract ions 
      trait_FA <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data_refined)) %>% 
        dplyr::filter(abs(mz-targ_df[i,"Mass_FA"])<(targ_df[i,"Mass_FA"]*ppm_error), abs(rt-targ_df[i, "RT"])<rt_error) 
      trait_H <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data_refined)) %>% 
        dplyr::filter(abs(mz-targ_df[i,"Mass_H"])<(targ_df[i,"Mass_H"]*ppm_error), abs(rt-targ_df[i, "RT"])<rt_error) 
      
      # Plot the chromatograms (FA)
      if(nrow(trait_FA) != 0){
        my_title <- paste(targ_df[i,"Cpd"], "[M+HCOO]-")
        peak_plot(res_PeakID$Peak_data_refined, trait_FA, my_title, nCore)
      }
      
      # Plot the chromatograms (H)
      if(nrow(trait_H) != 0){
        my_title <- paste(targ_df[i,"Cpd"], "[M-H]-")
        peak_plot(res_PeakID$Peak_data_refined, trait_H, my_title, nCore)
      }
      graphics.off()
    })  
  }
}

if(Align){
  # Peak identification check: Random set of peaks
  peakID_rand <- chromPeaks(res_PeakID$Peak_data_align) %>% 
    rownames() %>% 
    sample(., nrandom)
  peak_rand_df <- chromPeaks(res_PeakID$Peak_data_align)[peakID_rand,]
  
  if(untarg_plot){
    message("Plot random peaks (aligned)")
    # Plotting
    pdf(file.path(out_dir, "Plots", "XCMS_aligned_untargeted_peaks.pdf"), width = 12, height = 12)
    lapply(1:nrow(peak_rand_df), function(i){
      peak_ID <- rownames(peak_rand_df)[[i]]
      trait <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data_align)) %>% 
        dplyr::filter(abs(mz-peak_rand_df[i,"mz"])<(peak_rand_df[i,"mz"]*ppm_error), abs(rt-peak_rand_df[i, "rt"])<rt_error) 
      if(nrow(trait) != 0){
        peak_plot(res_PeakID$Peak_data_align, trait, peak_ID, nCore)
      }
    })
    graphics.off()
  }
  
  # Peak identification check: Targeted peaks
  if(targ_plot){
    message("Plot targeted peaks (aligned)")
    # Plotting
    dir.create(file.path(out_dir, "Plots", "Alignment_targeted_peaks"), showWarnings = FALSE)
    lapply(1:nrow(targ_df), function(i){
      # Prepare pdf
      trait_name <- targ_df[i,"Cpd"]
      fname <- file.path(out_dir, "Plots", "Alignment_targeted_peaks", paste0("Aligned_", i, "_", trait_name, ".png"))
      png(fname, width = 10, height = 4, units = 'in', res = 600)
      par(mfrow = c(1, 2))
      
      # Extract ions
      trait_FA <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data_align)) %>% 
        dplyr::filter(abs(mz-targ_df[i,"Mass_FA"])<(targ_df[i,"Mass_FA"]*ppm_error), abs(rt-targ_df[i, "RT"])<rt_error)
      trait_H <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data_align)) %>% 
        dplyr::filter(abs(mz-targ_df[i,"Mass_H"])<(targ_df[i,"Mass_H"]*ppm_error), abs(rt-targ_df[i, "RT"])<rt_error) 

      # Plot the chromatograms (FA)
      if(nrow(trait_FA) != 0){
        my_title <- paste(targ_df[i,"Cpd"], "[M+HCOO]-")
        peak_plot(res_PeakID$Peak_data_align, trait_FA, my_title, nCore)
      }
      # Plot the chromatograms (H)
      if(nrow(trait_H) != 0){
        my_title <- paste(targ_df[i,"Cpd"], "[M-H]-")
        peak_plot(res_PeakID$Peak_data_align, trait_H, my_title, nCore)
      }
      graphics.off()
    })
  }
  
  if(Align && ReFine){
    # Peak identification check: Random set of peaks
    peakID_rand <- chromPeaks(res_PeakID$Peak_data_refined_align) %>% 
      rownames() %>% 
      sample(., nrandom)
    peak_rand_df <- chromPeaks(res_PeakID$Peak_data_refined_align)[peakID_rand,]
    
    if(untarg_plot){
      message("Plot random peaks (refined + aligned)")
      
     # Plotting 
      pdf(file.path(out_dir, "Plots", "XCMS_refined_aligned_untargeted_peaks.pdf"), width = 12, height = 12)
      lapply(1:nrow(peak_rand_df), function(i){
        peak_ID <- rownames(peak_rand_df)[[i]]
        trait <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data_refined_align)) %>% 
          dplyr::filter(abs(mz-peak_rand_df[i,"mz"])<(peak_rand_df[i,"mz"]*ppm_error), abs(rt-peak_rand_df[i, "rt"])<rt_error) 
        if(nrow(trait) != 0){
          peak_plot(res_PeakID$Peak_data_refined_align, trait, peak_ID, nCore)
        }
      })
      graphics.off()
    }
    
    # Peak identification check: Targeted peaks
    if(targ_plot){
      message("Plot targeted peaks (refined + aligned)")
      
      # Create dir
      dir.create(file.path(out_dir, "Plots", "Refined_aligned_targeted_peaks"), showWarnings = FALSE)
      
      # Plotting 
      lapply(1:nrow(targ_df), function(i){
        
        # Prepare pdf 
        trait_name <- targ_df[i,"Cpd"]
        fname <- file.path(out_dir, "Plots", "Refined_aligned_targeted_peaks", paste0("Aligned_refined_", i, "_", trait_name, ".png"))
        png(fname, width = 10, height = 4, units = 'in', res = 600)
        par(mfrow = c(1, 2))
        
        # Extract ions
        trait_FA <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data_refined_align)) %>% 
          dplyr::filter(abs(mz-targ_df[i,"Mass_FA"])<(targ_df[i,"Mass_FA"]*ppm_error), abs(rt-targ_df[i, "RT"])<rt_error) 
        trait_H <- data.frame(xcms::chromPeaks(res_PeakID$Peak_data_refined_align)) %>% 
          dplyr::filter(abs(mz-targ_df[i,"Mass_H"])<(targ_df[i,"Mass_H"]*ppm_error), abs(rt-targ_df[i, "RT"])<rt_error) 
        
        # Plot the chromatograms (FA)
        if(nrow(trait_FA) != 0){
          my_title <- paste(targ_df[i,"Cpd"], "[M+HCOO]-")
          peak_plot(res_PeakID$Peak_data_refined_align, trait_FA, my_title, nCore)
        }
        # Plot the chromatograms (H)
        if(nrow(trait_H) != 0){
          my_title <- paste(targ_df[i,"Cpd"], "[M-H]-")
          peak_plot(res_PeakID$Peak_data_refined_align, trait_H, my_title, nCore)
        }
        graphics.off()
      })
    }
  }
}

# Save session info
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "Plots", paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")
tictoc::toc()