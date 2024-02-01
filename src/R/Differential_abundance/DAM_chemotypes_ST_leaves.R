#https://rawgit.com/andreasmock/MetaboDiff/master/vignettes/MetaboDiff_tutorial.html
suppressPackageStartupMessages({
  library(MetaboDiff)
  library(purrr)
  library(tidyverse)
  source("http://peterhaschke.com/Code/multiplot.R")
  source("~/Git/exAtlas/src/R/Differential_abundance/DAM_functions.R")
})

# Load data 
out_dir <- "~/Sara/exAtlas/res/Preprocessing_res/MS/DEM/Chemotypes_tissues"
untarg_int <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification/exAtlas_MS_untargeted_peak_intensities_tissue_normalized_singletons_rm.tsv")
exAtlas_LCMS_RNA_meta_info <- read.delim("~/Sara/exAtlas/res/RNAseq/analysis/DE/exAtlas_LCMS_RNA_meta_info.tsv") %>% 
  mutate(Chemotype = ifelse(Genotype %in% c(14, 34), "CN", NA),
         Chemotype = ifelse(Genotype %in% c(47, 60, 56), "AC", Chemotype),
         Chemotype = ifelse(Genotype %in% c(45, 41, 51, 76), "TL", Chemotype),
         Chemotype = as.factor(Chemotype),
         Genotype = as.factor(Genotype),
         Tissue_exp = as.factor(Tissue_exp))
write_tsv(exAtlas_LCMS_RNA_meta_info, file.path(out_dir, "exAtlas_LCMS_RNA_meta_info_tissue.tsv"))

# Set parameters
padj = 0.001
lfc = 0.5
my_tissues <- c("Shoot_tip_GH", "Leaf_D1_GH", "Leaf_D2_GH", "Leaf_D3_GH")

lapply(my_tissues, function(x){
  exAtlas_info <- exAtlas_LCMS_RNA_meta_info %>% filter(Tissue_exp %in% x) 
  sub_dir <- file.path(out_dir, x)
  dir.create(sub_dir, showWarnings = FALSE)
  
  # Prep tables
  untarg_int2 <- left_join(exAtlas_info %>% select(c(LCMS_sample, Chemotype)), untarg_int, by = c(LCMS_sample = "Sample"))
  untarg_int_dat <- untarg_int2 %>% dplyr::select(-c(Chemotype:Weight_mg)) %>% column_to_rownames("LCMS_sample")
  
  sample_meta <- untarg_int2 %>% 
    dplyr::select(c(LCMS_sample:Weight_mg)) %>% 
    column_to_rownames("LCMS_sample") %>% 
    mutate(Genotype = paste0("g", Genotype))
  
  meta_info <- data.frame(id = colnames(untarg_int_dat), 
                          RC = colnames(untarg_int_dat)) %>% column_to_rownames("id")
  
  (met <- create_mae(untarg_int_dat %>% t(), meta_info, sample_meta))
  (met = knn_impute(met,cutoff = 0.4))
  
  # Normalization
  (met <- normalize_met(met))
  
  # Quality control of normalization
  pdf(file.path(sub_dir, "Boxplot_QA_normalization.pdf"), width = 15, height = 8)
  quality_plot(met,
               group_factor = "Genotype",
               label_colors = c("darkseagreen","dodgerblue", "red", "blue", 
                                "green", "orange", "pink", "purple", "turquoise"))
  graphics.off()
  
  pdf(file.path(sub_dir, "PCA_tsne_genotype.pdf"), width = 15, height = 10)
  multiplot(
    pca_plot(met,
             group_factor = "Genotype",
             label_colors = c("darkseagreen","dodgerblue", "red", "blue", 
                              "green", "orange", "pink", "purple", "turquoise")),
    tsne_plot(met,
              group_factor = "Genotype",
              label_colors = c("darkseagreen","dodgerblue", "red", "blue", 
                               "green", "orange", "pink", "purple", "turquoise")),
    cols = 2)
  graphics.off()
  
  pdf(file.path(sub_dir, "PCA_tsne_chemotype.pdf"), width = 15, height = 10)
  multiplot(
    pca_plot(met,
             group_factor = "Chemotype",
             label_colors = c("darkseagreen","dodgerblue", "red")),
    tsne_plot(met,
              group_factor = "Chemotype",
              label_colors = c("darkseagreen","dodgerblue", "red")),
    cols = 2)
  graphics.off()
  
  # Differential abundance testing - tissue t-test comparison
  all_chem <- met@colData$Chemotype %>% unique()
  
  DEM_ttest <- lapply(all_chem, function(i){
    my_res <- lapply(all_chem, function(x) {
      if(i != x){
        diff_ttest(met, "Chemotype", c(x, i)) # using this test to make it comparable with our Differenital expression analysis (we want pairwise comparisons only)
      }
    })
    my_res %>% unlist(., recursive = FALSE)
  }) %>% unlist(., recursive = FALSE)
  
  # To match the DESeq logic (contrasts), we need to change the up and down regulation direction 
  DEM_lists <- lapply(DEM_ttest, function(x){
    my_res <- x
    sel <- my_res$adj_pval <= padj & abs(my_res$dm) >= lfc & ! is.na(my_res$adj_pval)
    RC_DEM <- list(all = my_res[sel,]$metabolite,
                   dn = my_res[sel & my_res$dm > 0,]$metabolite, # a positive dm means that it's up regulated in Y and down regulated in X (X_vs_Y). See MetaboDiff tutorial. 
                   up = my_res[sel & my_res$dm < 0,]$metabolite) 
    RC_DEM
  })
  names(DEM_lists) <- names(DEM_ttest)
  
  # Since it's duplicated, I will only keep half of the tests
  DEM_lists <- DEM_lists[1]
  
  # I also want TL to be the background (to match DESeq output)
  DEM_regulation <- list()
  for(x in names(DEM_lists)){
    my_contrast <- x 
    if(grepl("TL", my_contrast) & !endsWith(my_contrast, "TL")){
      if(grepl("CN", my_contrast)){
        new_contrast <- "ttest_Chemotype_CN_vs_TL"
      }
      if(grepl("AC", my_contrast)){
        new_contrast <- "ttest_Chemotype_AC_vs_TL"
      }
      DEM_regulation[[new_contrast]] <- list(
        all = DEM_lists[[my_contrast]]$all,
        dn = DEM_lists[[my_contrast]]$up,
        up = DEM_lists[[my_contrast]]$dn)
    } else {
      DEM_regulation[[my_contrast]] <- DEM_lists[[my_contrast]]
    }
  }
  
  # Save data
  save(met, DEM_ttest, DEM_regulation, file = file.path(sub_dir, paste0(paste0("exAtlas_DEM_chemotype_padj", padj,"_lfc", lfc,".RData"))))
})
