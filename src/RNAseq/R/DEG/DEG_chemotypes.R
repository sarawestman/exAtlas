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

pal = brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

# Load data
out_dir <- "~/Sara/exAtlas/res/RNAseq/analysis/DE/Chemotype_tissue"
my_tissues <- c("Shoot_tip_GH", "Leaf_D1_GH", "Leaf_D2_GH", "Leaf_D3_GH")

# Set parameters 
resOnly = TRUE
pval_thr = 0.001
lfc_thr = 0.5
filtering_type = "mean"

lapply(my_tissues, function(x){
  
  out_dir <- file.path(out_dir, x)
  dir.create(out_dir, showWarnings = FALSE)
  
  if(!resOnly){
    samples <- read_tsv("~/Sara/exAtlas/res/RNAseq/meta/RNA_seq_salmon_info.tsv", show_col_types = FALSE) %>%
      filter(Tissue_exp %in% x) %>% 
      mutate(Chemotype = ifelse(Genotype %in% c(14, 34), "CN", NA),
             Chemotype = ifelse(Genotype %in% c(47, 60, 56), "AC", Chemotype),
             Chemotype = ifelse(Genotype %in% c(45, 41, 51, 76), "TL", Chemotype),
             Chemotype = as.factor(Chemotype),
             Genotype = as.factor(Genotype),
             Tissue_exp = as.factor(Tissue_exp))
    write_tsv(samples, file.path(out_dir, "RNA_seq_sample_meta_tissue.tsv"))
    
    tx2gene <- suppressMessages(read_tsv("/mnt/picea/storage/reference/Populus-tremula/v2.2/annotation/tx2gene.tsv.gz",
                                         col_types=cols(TXID = col_character(), GENEID = col_character())))
    
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
      design = ~ Chemotype) 
    colnames(dds) <- samples$sample
    save(dds, file = file.path(out_dir, "dds.rda"))
    
    # vst 
    vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
    vst <- assay(vsd)
    vst <- vst - min(vst)
    save(vst, file = file.path(out_dir, "vst-aware.rda"))
    write_delim(as.data.frame(vst) %>% rownames_to_column("ID"), file.path(out_dir," vst-aware.tsv"))
  } else{
    load(file.path(out_dir, "dds.rda"))
    load(file.path(out_dir, "vst-aware.rda"))
  }
  
  dds$Chemotype <- relevel(dds$Chemotype, "TL")
  dds <- DESeq(dds)
  resultsNames(dds)
  
  # Plot dispersion
  pdf(file.path(out_dir, "Dispersion_plot.pdf"))
  plotDispEsts(dds)
  graphics.off()
  
  # Make output dir 
  if(filtering_type == "mean"){
    filtering_sel <- "NULL"
    out_dir <- file.path(out_dir, "mean")
    dir.create(out_dir, showWarnings = FALSE)
  } else if (filtering_type == "median"){
    filtering_sel = "median"
    out_dir <- file.path(out_dir, "median")
    dir.create(out_dir, showWarnings = FALSE)
  }
  
  sub_dir <- file.path(out_dir, paste0("padj", pval_thr, "_lfc", lfc_thr))
  dir.create(sub_dir, showWarnings = FALSE)
  
  # Extract data 
  lapply(resultsNames(dds)[2:3], function(i){
    contrasts_sel <- i
    DE_res <- extract_results(dds, 
                              vst, 
                              padj = pval_thr,
                              lfc = lfc_thr,
                              contrasts_sel,
                              filter = filtering_sel,
                              default_dir = sub_dir,
                              default_prefix = contrasts_sel,
                              labels = dds$Genotype,
                              plot = TRUE, 
                              verbose = TRUE)
    DE_res$Contrast <- contrasts_sel
    save(DE_res, file = file.path(sub_dir, paste0(contrasts_sel,".rda")))
  })
})
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
