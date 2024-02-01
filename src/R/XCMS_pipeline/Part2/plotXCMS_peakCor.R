#######################################
# # Plot peak correspondance
#######################################
source("~/Git/exAtlas/src/R/XCMS_pipeline/XCMS_functions.R")
suppressPackageStartupMessages({
  library(xcms)
  library(randomcoloR)
  library(RColorBrewer)
  library(magrittr)
  library(SummarizedExperiment)
})
tictoc::tic()

################################################           
# Load files  
################################################
# Output directory
out_dir <- "~/Sara/exAtlas/res/XCMS/MS/Peak_correspondance"
in_file <- "~/Sara/exAtlas/res/XCMS/MS/Peak_correspondance/exAtlas_MS_XCMS_peak_correspondance.RData"
targ_file <- "~/Sara/exAtlas/meta/Targeted_metabolites/SPG_targeted_complete_list_adducts.tsv"
  
# Output file prefix 
file_name <- "exAtlas_MS"

# Table with targeted compounds 
targ_df <- read.delim(targ_file)

# XCMS PeakCor object 
load(in_file)

################################################           
# Set plotting parameters   
################################################
# Which plots would you like to do?
targ_plot = TRUE
untarg_plot = TRUE
PCA_plot = TRUE

# XCMS object to plot 
XCMS_obj_sel = "PeakID_obj_gr_filled" # use PeakID_obj_gr in the feature (inc. meta data) - OR readr was all I needed? 

# Column for ChromPeak density plot coloring
samp_coloring = "Sample"

# PCA plot coloring 
PCA_coloring = "Tissue_exp"
PCA_lab = "Genotype_nr"

# Set ppm/rt error for targeted peaks 
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
  message(paste("Output file prefix:", file_name))
}

# Set parallelisation - nr of cores
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1) {
  stop("Nr of cores for parallelisation must be supplied", call.=FALSE)
} else {
  nCore <- args[1]
  message("Nr of cores: ", nCore)
}

################################################           
# Plot peaks  
################################################
# Select XCMS objects for plotting 
message(paste("XCMS object selected for plotting:", XCMS_obj_sel))
XCMS_obj <- res_cor[[XCMS_obj_sel]]

# Set colors for plotting
nb.cols <- length(pData(XCMS_obj)[[samp_coloring]] %>% unique())
sample_colors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
names(sample_colors) <- pData(XCMS_obj)[[samp_coloring]] %>% unique()

# Plot targeted peaks 
if(targ_plot){
  message("Plot targeted chromatogram peaks")
  pdf(file.path(out_dir, "Plots", paste0(file_name, "_ChromPeak_density_plots_targeted_peaks.pdf")), width = 12, height = 12)
  lapply(1:nrow(targ_df), function(i){
    trait_name <- targ_df[i,"Cpd"]
    trait_H <- data.frame(xcms::chromPeaks(XCMS_obj)) %>% 
      dplyr::filter(abs(mz-targ_df[i,"Mass_H"])<(targ_df[i,"Mass_H"]*ppm_error),abs(rt-targ_df[i, "RT"])<rt_error) 
    trait_FA <- data.frame(xcms::chromPeaks(XCMS_obj)) %>% 
      dplyr::filter(abs(mz-targ_df[i,"Mass_FA"])<(targ_df[i,"Mass_FA"]*ppm_error),abs(rt-targ_df[i, "RT"])<rt_error) 
  
    # H loss
    if(nrow(trait_H) != 0){
      my_title <- paste(trait_name, "[M-H]-")
      peak_plot(XCMS_obj, trait_H, my_title)
      peak_ChromPeakDensity_plot(XCMS_obj, trait_H, sample_colors, res_cor$PeakDensity_Param, my_title, nCore)
    }
    # FA adduct
    if(nrow(trait_FA) != 0){
      my_title <- paste(trait_name, "[M+COOH]-")
      peak_plot(XCMS_obj, trait_FA, my_title)
      peak_ChromPeakDensity_plot(XCMS_obj, trait_FA, sample_colors, res_cor$PeakDensity_Param, my_title, nCore)
    }
    })
  graphics.off()
  }

# Plot random features 
if(untarg_plot){
  message("Plot random peaks")
  feature_rand <- sample(featureValues(XCMS_obj) %>% rownames(), 100)
  pdf(file.path(out_dir, "Plots", paste0(file_name, "_Feature_chromatogram_random_peaks.pdf")), width = 12, height = 12)
  lapply(feature_rand, function(i){
    feature_chroms <- featureChromatograms(XCMS_obj, 
                                           features = i)
    plot(feature_chroms, 
         col = sample_colors[chromPeaks(feature_chroms)[, "sample"]],
         peakBg = sample_colors[chromPeaks(feature_chroms)[, "sample"]])
    })
  graphics.off()
}

################################################           
# Plot PCA  
################################################
## Extract the features and log2 transform them
if(PCA_plot){
  message("Plot PCA")
  ft_ints <- log2(assay(res_cor$Peak_overview, "raw_filled"))

  ## Perform the PCA omitting all features with an NA in any of the samples. Also, the intensities are mean centered.
  pc <- prcomp(t(na.omit(ft_ints)), center = TRUE)
  pcSummary <- summary(pc)

  # Set colors for plotting 
  group_colors <- randomColor(nb.cols, luminosity = "random")
  names(group_colors) <- pData(XCMS_obj)[[PCA_coloring]] %>% unique()
  cols <- group_colors[pData(XCMS_obj)[[PCA_coloring]]]

  # PLot PCA 
  pdf(file.path(out_dir, "Plots", paste0(file_name, "_PCA.pdf")))
  plot(pc$x[, 1], pc$x[,2], pch = 21, main = "",
       xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
                                     digits = 3), " % variance"),
       ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
                                     digits = 3), " % variance"),
      col = "darkgrey", bg = cols, cex = 2)
  grid()
  text(pc$x[, 1], pc$x[,2], labels = pData(XCMS_obj)[[PCA_lab]], col = "darkgrey",
       pos = 3, cex = 1)
  graphics.off()
}

# Save session info
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "Plots", paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")
tictoc::toc()