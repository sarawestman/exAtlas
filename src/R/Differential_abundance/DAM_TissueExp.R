suppressPackageStartupMessages({
  library(MetaboDiff)
  library(purrr)
  library(tidyverse)
  source("http://peterhaschke.com/Code/multiplot.R")
  source("~/Git/exAtlas/src/R/Differential_abundance/extra_func.R")
})

# Load data 
out_dir <- "~/Sara/exAtlas/res/Preprocessing_res/MS/DEM/Tissue_exp"
untarg_int <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification/exAtlas_MS_untargeted_peak_intensities_tissue_normalized_singletons_rm.tsv")
exAtlas_LCMS_RNA_meta_info <- read.delim("~/Sara/exAtlas/res/RNAseq/analysis/DE/exAtlas_LCMS_RNA_meta_info.tsv")

# Prep tables
untarg_int2 <- left_join(exAtlas_LCMS_RNA_meta_info %>% select(c(LCMS_sample)), untarg_int, by = c(LCMS_sample = "Sample"))
untarg_int_dat <- untarg_int2 %>% dplyr::select(-c(Injection:Weight_mg)) %>% column_to_rownames("LCMS_sample")

sample_meta <- untarg_int2 %>% 
  dplyr::select(c(LCMS_sample:Weight_mg)) %>% 
  column_to_rownames("LCMS_sample") %>% 
  mutate(Genotype = paste0("g", Genotype))

meta_info <- data.frame(id = colnames(untarg_int_dat), 
                        RC = colnames(untarg_int_dat)) %>% column_to_rownames("id")

(met <- create_mae(untarg_int_dat %>% t(), meta_info, sample_meta))
(met = knn_impute(met, cutoff = 0.4))

# Normalization
(met <- normalize_met(met))

LCMS_norm <- met@ExperimentList$norm_imputed %>% assay() %>% t() %>% as.data.frame() %>% rownames_to_column("Sample")
write_tsv(LCMS_norm, file.path(out_dir, "exAtlas_LCMS_vsn_normalized_intensity.tsv"))

# Quality control of normalization
pdf(file.path(out_dir, "Boxplot_QA_normalization.pdf"), width = 15, height = 8)
quality_plot(met,
             group_factor = "Tissue_exp",
             label_colors = c("darkseagreen","dodgerblue", "red", 
                              "blue", "green", "orange", 
                              "pink", "purple", "darkkhaki"))
graphics.off()

# Unsupervised analysis
pdf(file.path(out_dir, "PCA_tsne.pdf"), width = 15, height = 10)
multiplot(
  pca_plot(met,
           group_factor = "Tissue_exp",
           label_colors = c("darkseagreen","dodgerblue", "red", 
                            "blue", "green", "orange", 
                            "pink", "purple", "darkkhaki")),
  tsne_plot(met,
            group_factor = "Tissue_exp",
            label_colors = c("darkseagreen","dodgerblue", "red", 
                             "blue", "green", "orange", 
                             "pink", "purple", "darkkhaki")),
  cols=2)
graphics.off()

pdf(file.path(out_dir, "PCA_tsne2.pdf"), width = 15, height = 10)
multiplot(
  pca_plot(met,
           group_factor = "Genotype",
           label_colors = c("darkseagreen","dodgerblue", "red", "blue", "green", "orange",
                          "black", "pink", "purple")),
  tsne_plot(met,
            group_factor="Genotype",
            label_colors = c("darkseagreen","dodgerblue", "red", "blue", "green", "orange",
                           "black", "pink", "purple")),
  cols=2)
graphics.off()

# Differential abundance testing - tissue t-test comparison
all_tissues <- met@colData$Tissue_exp %>% unique()

DEM_ttest <- lapply(all_tissues, function(i){
  my_res <- lapply(all_tissues, function(x) {
    if(i != x){
      diff_test2(met, "Tissue_exp", c(x, i))
    }
  })
  my_res %>% unlist(., recursive = FALSE)
}) %>% unlist(., recursive = FALSE)

save(met, DEM_ttest, file = file.path(out_dir, "exAtlas_DEM_tissue_exp.RData"))

#volcano_plot2(DEM_ttest, 
#              comp = "ttest_Tissue_exp_Leaf_D3_GH_vs_Leaf_D2_GH",
#              label_colors = c("dodgerblue","firebrick4"),
#              dm_cutoff = 0.5,
#              p_adjust = TRUE)


