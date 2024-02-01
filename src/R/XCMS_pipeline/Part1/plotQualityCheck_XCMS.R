#########################################################
# XCMS - Part 1: Plot XCMS peak identification results
#########################################################
source("~/Git/exAtlas/src/R/XCMS_pipeline/XCMS_functions.R")
suppressPackageStartupMessages({
  library(xcms)
  library(RColorBrewer)
  library(pander)
  library(magrittr)
  library(pheatmap)
  library(SummarizedExperiment)
  library(BiocParallel)
})
tictoc::tic()

################################################           
# Load files  
################################################
# Output directory
out_dir <- "~/Sara/exAtlas/res/XCMS/MS/Peak_identification"
in_file <- "~/Sara/exAtlas/res/XCMS/MS/Peak_identification/exAtlas_MS_XCMS_peak_identification.RData"
targ_file <- "~/Sara/exAtlas/meta/Targeted_metabolites/SPG_targeted_complete_list_adducts.tsv"

# Set prefix 
file_name <- "exAtlas_MS"

# Targeted data 
targ_df <- read.delim(targ_file, header = TRUE)

# Load xcms data 
load(in_file)

################################################           
# Set XCMS parameters  
################################################
# Plot results from intensity filtering? 
ReFine = TRUE

# Plot results from alignment? 
Align = FALSE
sampleGroupsType = res_PeakID$Meta_info$Type
sampleGroupsType_sel = "QC"

# Quality check plots: coloring groups
gr_coloring = "Tissue_exp" 
gr_coloring_spec = paste0(res_PeakID$Meta_info$Genotype_nr, "_Inj", res_PeakID$Meta_info$Injection) # needs to be unique labels! 

######################################
# Check assumptions
######################################
# Check if gr_coloring_spec contains unique labels  
if(!length(unique(gr_coloring_spec)) == length(gr_coloring_spec)){
  stop("'gr_coloring_spec' must contain unique lables")
}

# Check if output directory exists 
if(!dir.exists(out_dir)){
  stop(paste("Output directory does not exist:", out_dir))
} else{
  message(paste("Input file:", in_file))
  message(paste("Targeted file used:", targ_file))
  dir.create(file.path(out_dir, "Quality_check"), showWarnings = FALSE)
  message(paste("Output directory:", file.path(out_dir, "Quality_check")))
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

######################################
# Plot - Quality checks
######################################
# Plot the tic
message("Plot TIC")
chr <- xcms::chromatogram(res_PeakID$Raw_data, 
                    aggregationFun = "max",
                    BPPARAM = BiocParallel::MulticoreParam(workers = nCore))
pdf(file.path(out_dir, paste0("TIC_", file_name,".pdf")), width = 12)
plot(chr)
graphics.off()

# Plot alignment tic 
if(Align){
  # Set colors
  message("Plot alignment")
  res_PeakID$Peak_data_align$Type <- sampleGroupsType
  clrs <- rep("#00000040", length(res_PeakID$Peak_data_align$Type))
  clrs[res_PeakID$Peak_data_align$Type == sampleGroupsType_sel] <- c("#00ce0080")
  
  # Create plot
  pdf(file.path(out_dir, "Quality_check", paste0("TIC_alignment_", file_name,".pdf")), width = 12)
  par(mfrow = c(2, 1), mar = c(4, 4.5, 1, 0.5))
  plot(xcms::chromatogram(res_PeakID$Peak_data_align, 
                    aggregationFun = "sum",
                    BPPARAM = BiocParallel::MulticoreParam(workers = nCore)),
       col = clrs, peakType = "none")
  plotAdjustedRtime(res_PeakID$Peak_data_align, col = clrs, peakGroupsPch = 1,
                    peakGroupsCol = "#00ce0040")
  graphics.off()
}
if(Align && ReFine){
  # Set colors
  message("Plot alignment (refined + aligned)")
  res_PeakID$Peak_data_refined_align$Type <- sampleGroupsType
  clrs <- rep("#00000040", length(res_PeakID$Peak_data_refined_align$Type))
  clrs[res_PeakID$Peak_data_refined_align$Type == sampleGroupsType_sel] <- c("#00ce0080")
  
  # Create plot
  pdf(file.path(out_dir, "Quality_check", paste0("TIC_alignment_refined_", file_name,".pdf")), width = 12)
  par(mfrow = c(2, 1), mar = c(4, 4.5, 1, 0.5))
  plot(xcms::chromatogram(res_PeakID$Peak_data_refined_align, 
                    aggregationFun = "sum", 
                    BPPARAM = BiocParallel::MulticoreParam(workers = nCore)),
       col = clrs, peakType = "none")
  plotAdjustedRtime(res_PeakID$Peak_data_refined_align, 
                    col = clrs, 
                    peakGroupsPch = 1,
                    peakGroupsCol = "#00ce0040")
  graphics.off()
}

# Below we create boxplots representing the distribution of total ion currents per file. 
# Such plots can be very useful to spot problematic or failing MS runs.
tc <- split(tic(res_PeakID$Raw_data), f = fromFile(res_PeakID$Raw_data))
nb.cols <- length(pData(res_PeakID$Raw_data)[[gr_coloring]] %>% unique())
group_colors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
names(group_colors) <- pData(res_PeakID$Raw_data)[[gr_coloring]] %>% unique()

# Plot boxplots 
message("Plot boxplot")
pdf(file.path(out_dir, "Quality_check", paste0("Boxplot_TIC_", file_name, ".pdf")), width = 30)
par(mar = c(10, 5, 5, 5))
boxplot(tc, 
        col =  group_colors[pData(res_PeakID$Raw_data)[[gr_coloring]]],
        ylab = "intensity", 
        names = pData(res_PeakID$Raw_data)[[gr_coloring]],
        las = 2,
        main = "Total ion current") 
graphics.off()

# Quality check: pheatmap to identify outliers 
chr_bin <- MSnbase::bin(chr, binSize = 2)

# Calculate correlation on the log2 transformed base peak intensities and define which phenodata columns should be highlighted in the plot
cormat <- cor(log2(do.call(cbind, lapply(chr_bin, intensity))))
ann <- data.frame(group = pData(res_PeakID$Raw_data)[[gr_coloring]])
rownames(ann) <- colnames(cormat) <- rownames(cormat) <- gr_coloring_spec
start_colnames <- colnames(cormat)

# Check for NA in correlation matrix 
if(any(is.na(cormat))){
  cormat <- cormat[-which(is.na(cormat)),-which(is.na(cormat))]
  new_colnames <- colnames(cormat)
  message(paste("Sample(s) removed from pheatmap due to NA:", start_colnames[!start_colnames %in% new_colnames]))
  message("Saving pheatmap data for quality check")
  save(chr, chr_bin, cormat, ann, group_colors, file = file.path(out_dir, "Pheatmap_data_control_NA.RData"))
}

# Plot pheatmap 
message("Plot pheatmap")
pdf(file.path(out_dir, "Quality_check", paste0("Pheatmap_", file_name, ".pdf")), width = 30, height = 30)
pheatmap(cormat, 
         border_color = NA,
         annotation = ann,
         annotation_color = list(group = group_colors))
graphics.off()

# Save session info
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "Quality_check", paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")
tictoc::toc()