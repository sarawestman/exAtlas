suppressPackageStartupMessages({
  library(reshape2)
  library(rrvgo)
  library(topGO) #https://gist.github.com/slavailn/dcff753cc32fd9e053590308c831b057
  library(plyr)
  library(tidyverse)
})

# Load data  
out_dir <- "~/Sara/exAtlas/res/RNAseq/analysis/DE/Tissue_vs_everything_else"
DE_files <- list.files("~/Sara/exAtlas/res/RNAseq/analysis/DE/Tissue_vs_everything_else", "_Yes_vs_No.rda",  recursive = TRUE, full.names = TRUE)[1:7]
geneID2GO <- readMappings("/mnt/picea/projects/aspseq/nstreet/swasp/Sara/Genome_analyses/Annotations/potra2go.tsv")
geneNames <- names(geneID2GO)

mode <- "dn"
DEG_res_sel <- lapply(DE_files, function(x){
  load(x)
  DE_res[[mode]]
})
names(DEG_res_sel)  <- gsub("(?<!^)([[:upper:]])", " \\1", str_remove_all(basename(dirname(DE_files)), paste(c("is", "GH"), collapse = "|")), perl = TRUE)

# Select data 
DEG_sel <- DEG_res_sel[c(6,2:4,7,1,5)]

# Run enrichment 
enr_res <- lapply(names(DEG_sel), function(x){
  message(x)
  ItemsList_genes <- venn::venn(DEG_sel, show.plot = FALSE)
  ItemsList_sel_genes <- attributes(ItemsList_genes)$intersections[x][[1]]
  
  geneList <- factor(as.integer(geneNames %in% ItemsList_sel_genes))
  names(geneList) <- geneNames
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  allGO <- usedGO(GOdata)
  resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata, pval = resultFis, orderBy = 'pval', topNodes = length(allGO)) %>% 
    mutate(pval = as.numeric(pval),
           padj = round(p.adjust(pval, method = "BH"), digits = 4)) %>% 
    arrange(padj)
})
names(enr_res) <- names(DEG_sel)

# Plot 
treemap::treemap(
  enr_res[[1]] %>% filter(padj < 0.05),
  index = "Term", 
  vSize = "pval",
  vColor = "pval")

# Run revigo
enr_revigo <- lapply(names(enr_res), function(x){
  allResSig <- enr_res[[x]] %>% filter(padj < 0.05)
  if(nrow(allResSig) >= 2){
    simMatrix <- calculateSimMatrix(allResSig$GO.ID,
                                    orgdb = "org.At.tair.db",
                                    ont = "BP",
                                    method = "Rel")

    scores <- setNames(-log10(allResSig$pval), allResSig$GO.ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold = 0.9,
                                    orgdb = "org.At.tair.db")
    }
  })
names(enr_revigo) <- names(enr_res)

enr_revigo <- plyr::compact(enr_revigo)

# Save results 
save(enr_res, enr_revigo, file = file.path(out_dir, paste0("TopGO_revigo_res_GH_", mode,".RData")))

# Make dir
sub_dir <- file.path(out_dir, "Plots")
dir.create(sub_dir, showWarnings = FALSE)
sub_dir2 <- file.path(sub_dir, paste0("wordcloud_topGO_", mode))
dir.create(sub_dir2, showWarnings = FALSE)

# Produce plots 
lapply(names(enr_revigo), function(x){
  reducedTerms2 <- enr_revigo[[x]] %>% dplyr::select(c(parentTerm, score)) %>% group_by(parentTerm) %>% top_n(1, score) %>% distinct()
  pdf(file.path(sub_dir2, paste0("WC_", x, "_topGO_parentTerm.pdf")))
  wordcloud::wordcloud(words = reducedTerms2$parentTerm, freq = reducedTerms2$score, min.freq = 0,           
                       max.words = 200, random.order = FALSE, rot.per = 0, scale = c(1,.4),          
                       colors = RColorBrewer::brewer.pal(8, "Dark2"))
  graphics.off()
})

writeLines(capture.output(sessionInfo()), file.path(sub_dir2, "sessionInfo.txt"))

