suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(RColorBrewer)
  library(tximport)
  library(tidyverse)
  library(VennDiagram)
})
suppressMessages({
  #source("~/src/UPSCb-common/Rtoolbox/src/plotEnrichedTreemap.R")
  source("~/Git/exAtlas/src/RNAseq/R/DE_functions.R")
  source("~/src/UPSCb-common/src/R/featureSelection.R")
  source("~/src/UPSCb-common/src/R/volcanoPlot.R")
  source("~/src/UPSCb-common/src/R/gopher.R")
  source("~/Git/exAtlas/src/RNAseq/R/plotEnrichedTreemap.R")
})

# Colors 
pal = brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

# Load data
out_dir <- "~/Sara/exAtlas/res/RNAseq/analysis/DE/Tissue_vs_everything_else"

samples <- read_tsv("~/Sara/exAtlas/res/RNAseq/meta/RNA_seq_salmon_info.tsv", show_col_types = FALSE) %>% 
  filter(!Tissue_exp %in% c("Bark_Savar", "Buds_Savar")) %>% 
  mutate(Genotype = as.factor(Genotype), 
         Tissue_exp = as.factor(Tissue_exp), 
         isRoot = as.factor(ifelse(Tissue_exp == "Roots_GH", "Yes", "No")),
         isLeafD1 = as.factor(ifelse(Tissue_exp == "Leaf_D1_GH", "Yes", "No")),
         isLeafD2 = as.factor(ifelse(Tissue_exp == "Leaf_D2_GH", "Yes", "No")),
         isLeafD3 = as.factor(ifelse(Tissue_exp == "Leaf_D3_GH", "Yes", "No")),
         isShootTip = as.factor(ifelse(Tissue_exp == "Shoot_tip_GH", "Yes", "No")),
         isBudsGH = as.factor(ifelse(Tissue_exp == "Buds_GH", "Yes", "No")),
         isStem = as.factor(ifelse(Tissue_exp == "Stem_GH", "Yes", "No")))
write_tsv(samples, file.path(out_dir, "RNA_seq_sample_meta_tissue.tsv"))

tx2gene <- suppressMessages(read_tsv("/mnt/picea/storage/reference/Populus-tremula/v2.2/annotation/tx2gene.tsv.gz", col_types = cols(TXID = col_character(), GENEID = col_character())))

txi <- suppressMessages(tximport(files = samples$salmon_path,
                                 type = "salmon",
                                 tx2gene = tx2gene))

# Merge technical reps 
txi$counts <- sapply(split.data.frame(t(txi$counts),samples$Tube_nr),colSums)
txi$length <- sapply(split.data.frame(t(txi$length),samples$Tube_nr),colMaxs)

#' # Counts are now in alphabetic order, check and reorder if necessary
#stopifnot(colnames(txi$counts) == samples$Tube_nr)
samples <- samples[match(colnames(txi$counts),samples$Tube_nr),]
stopifnot(colnames(txi$counts) == samples$Tube_nr)

# On merged technical replicates? 
dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = samples,
  design = ~ Tissue_exp) 
colnames(dds) <- samples$sample
save(dds, file = file.path(out_dir, "dds.rda"))

# vst 
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
save(vst, file = file.path(out_dir, "vst-aware.rda"))
write_delim(as.data.frame(vst) %>% rownames_to_column("ID"), file.path(out_dir," vst-aware.tsv"))

# Differential expression analysis
my_baseline <- "No"
formulas_col <- c("isRoot","isLeafD1","isLeafD2","isLeafD3","isShootTip","isBudsGH","isStem")
filtering_type = "mean"
pval_thr = 0.001
lfc_thr = 0.5

lapply(formulas_col, function(x){
  # Run model 
  my_formula <- x
  dds[[x]] <- relevel(dds[[x]], my_baseline)
  design(dds) <- formula(paste0("~", my_formula))
  dds <- DESeq(dds)
  resultsNames(dds)
  
  # Make dir 
  DE_dir <- file.path(out_dir, my_formula)
  dir.create(file.path(DE_dir), showWarnings = FALSE)
  
  # Set filtering 
  if(filtering_type == "mean"){
    filtering_sel <- "NULL"
    out_dir <- file.path(DE_dir, "mean")
    dir.create(out_dir, showWarnings = FALSE)
  } else if (filtering_type == "median"){
    filtering_sel = "median"
    out_dir <- file.path(DE_dir, "median")
    dir.create(out_dir, showWarnings = FALSE)
  }
  
  sub_dir <- file.path(out_dir, paste0("padj", pval_thr, "_lfc", lfc_thr))
  dir.create(sub_dir, showWarnings = FALSE)
  
  # Plot dispersion
  pdf(file.path(DE_dir, "Dispersion_plot.pdf"))
  plotDispEsts(dds)
  graphics.off()
  
  # Extract data 
  contrasts_sel <- resultsNames(dds)[startsWith(resultsNames(dds), my_formula[[1]])]
  DE_res <- extract_results(dds, 
                            vst, 
                            padj = pval_thr,
                            lfc = lfc_thr,
                            contrasts_sel,
                            filter = filtering_sel,
                            default_dir = sub_dir,
                            default_prefix = contrasts_sel,
                            labels = dds$Tissue_exp_geno,
                            plot = TRUE, 
                            verbose = TRUE)
  DE_res$Contrast <- contrasts_sel
  
  save(DE_res, file = file.path(DE_dir, paste0(contrasts_sel,".rda")))
})
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))