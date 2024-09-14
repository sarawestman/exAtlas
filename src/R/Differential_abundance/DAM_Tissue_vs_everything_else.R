suppressPackageStartupMessages({
  library(MetaboDiff)
  library(purrr)
  library(tidyverse)
  source("http://peterhaschke.com/Code/multiplot.R")
})

# Load data 
out_dir <- "~/Sara/exAtlas/res/Preprocessing_res/MS/DEM/Tissue_vs_everything_else"
untarg_int <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification/exAtlas_MS_untargeted_peak_intensities_tissue_normalized_singletons_rm.tsv")
exAtlas_LCMS_RNA_meta_info <- read.delim("~/Sara/exAtlas/res/RNAseq/analysis/DE/exAtlas_LCMS_RNA_meta_info.tsv") %>% 
  filter(!Tissue_exp %in% c("Bark_Savar", "Buds_Savar")) %>% 
  select(-c(Chemotype, Tissue_chemotype)) %>% 
  mutate(isLeafD1 = as.factor(ifelse(Tissue_exp == "Leaf_D1_GH", "Yes", "No")),
         isLeafD2 = as.factor(ifelse(Tissue_exp == "Leaf_D2_GH", "Yes", "No")),
         isLeafD3 = as.factor(ifelse(Tissue_exp == "Leaf_D3_GH", "Yes", "No")),
         isShootTip = as.factor(ifelse(Tissue_exp == "Shoot_tip_GH", "Yes", "No")),
         isBudsGH = as.factor(ifelse(Tissue_exp == "Buds_GH", "Yes", "No")),
         isStem = as.factor(ifelse(Tissue_exp == "Stem_GH", "Yes", "No"))
  )
write_tsv(exAtlas_LCMS_RNA_meta_info, file.path(out_dir, "exAtlas_LCMS_RNA_meta_info_tissue.tsv"))

# Prep tables
untarg_int2 <- left_join(exAtlas_LCMS_RNA_meta_info %>% select(c(LCMS_sample, isRoot:isStem)), untarg_int, by = c(LCMS_sample = "Sample"))
untarg_int_dat <- untarg_int2 %>% dplyr::select(-c(isRoot:Weight_mg)) %>% column_to_rownames("LCMS_sample")

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

# Quality control of normalization
pdf(file.path(out_dir, "Boxplot_QA_normalization.pdf"), width = 15, height = 8)
quality_plot(met,
             group_factor = "Tissue_exp",
             label_colors = c("darkseagreen","dodgerblue", "red", 
                              "blue", "green", "orange", 
                              "pink"))
graphics.off()

# Unsupervised analysis
pdf(file.path(out_dir, "PCA_tsne_tissue.pdf"), width = 15, height = 10)
multiplot(
  pca_plot(met,
           group_factor="Tissue_exp",
           label_colors = c("darkseagreen","dodgerblue", "red", 
                            "blue", "green", "orange", 
                            "pink")),
  tsne_plot(met,
            group_factor="Tissue_exp",
            label_colors = c("darkseagreen","dodgerblue", "red", 
                             "blue", "green", "orange", 
                             "pink")),
  cols=2)
graphics.off()

pdf(file.path(out_dir, "PCA_tsne_genotype.pdf"), width = 15, height = 10)
multiplot(
  pca_plot(met,
           group_factor="Genotype",
           label_colors=c("darkseagreen","dodgerblue", "red", "blue", "green", "orange",
                          "black", "pink", "purple")),
  tsne_plot(met,
            group_factor="Genotype",
            label_colors=c("darkseagreen","dodgerblue", "red", "blue", "green", "orange",
                           "black", "pink", "purple")),
  cols=2)
graphics.off()

# Hypothesis testing
test_col <- c("isRoot", "isLeafD1", "isLeafD2", "isLeafD3", "isShootTip", "isBudsGH", "isStem")
DEM_test <- lapply(test_col, function(x){
  my_res <- diff_test(met, group_factors = c(x))
  metadata(my_res)[1]
})
names(DEM_test) <- test_col

# Filtering 
padj = 0.001
lfc = 0.5

DEM_lists <- lapply(DEM_test, function(x){
  my_res <- x[[1]]
  sel <-my_res$adj_pval <= padj & abs(my_res$dm) >= lfc & ! is.na(my_res$adj_pval)
  RC_DEM <- list(all = my_res[sel,]$metabolite,
                 dn = my_res[sel & my_res$dm > 0,]$metabolite, # a positive dm means that it's up regulated in Y and down regulated in X (X_vs_Y). See MetaboDiff tutorial. 
                 up = my_res[sel & my_res$dm < 0,]$metabolite) 
  RC_DEM
  })
names(DEM_lists) <- names(DEM_test)

# Save data
save(met, DEM_test, DEM_lists, file = file.path(out_dir, "exAtlas_LCMS_DEM_Tissue_vs_everything_else.RData"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
