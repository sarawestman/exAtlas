#' ---
#' title: "exAtlas Biological QA"
#' subtitle: "exAtlas"
#' author: "Sara"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    fig_width: 9
#'    fig_height: 6
#'    toc: true
#'    number_sections: true
#'    toc_depth: 3
#'    toc_float:
#'      collapsed: TRUE
#'      smooth_scroll: TRUE
#'    code_folding: hide
#'    theme: "flatly"
#'    highlight: pygments
#'    includes:
#'      before_body: header.html
#'      after_body: footer.html
#'    css: style.css
#' ---
#' 
#' <hr />
#' &nbsp;
#' 
#' # Setup
#' This section and the next are relevant for reproducibility purposes. For results, please skip to section 3 (Quality Control)
#' 
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(emoji)
  library(gplots)
  library(gtools)
  library(here)
  library(hyperSpec)
  library(limma)
  library(magrittr)
  library(parallel)
  library(patchwork)
  library(pheatmap)
  library(plotly)
  library(pvclust)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source("~/Git/exAtlas/src/RNAseq/R/featureSelection.R")
out_dir <- "~/Sara/exAtlas/res/RNAseq/analysis/salmon"
dir.create(file.path(out_dir, "QA_no_merged_technical_reps"),showWarnings=FALSE,recursive=TRUE)

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' # Data
#' * Sample information
#' ```{r Instructions1,eval=FALSE,echo=FALSE}
#' # The csv file should contain the sample information, including the sequencing file name, 
#' # any relevant identifier, and the metadata of importance to the study design
#' # as columns, e.g. the SamplingTime for a time series experiment
#'  ```
samples <- read_tsv("~/Sara/exAtlas/res/RNAseq/meta/RNA_seq_salmon_info.tsv")

#' * tx2gene translation table
#' ```{r Instructions2,eval=FALSE,echo=FALSE}
#' # This file is necessary if your species has more than one transcript per gene.
#' #
#' # It should then contain two columns, tab delimited, the first one with the transcript
#' # IDs and the second one the corresponding
#' #
#' # If your species has only one transcript per gene, e.g. Picea abies v1, then
#' # comment the next line
#' ```
tx2gene <- suppressMessages(read_tsv("/mnt/picea/storage/reference/Populus-tremula/v2.2/annotation/tx2gene.tsv.gz",
                                     col_types = cols(TXID = col_character(), 
                                                      GENEID = col_character())))
#' * Raw data
#filelist <- samples$salmon_path

#' * Sanity check to ensure that the data is sorted according to the sample info
#' ```{r Instructions3,eval=FALSE,echo=FALSE}
#' # This step is to validate that the salmon files are inthe same order as 
#' # described in the samples object. If not, then they need to be sorted
#' ````
#stopifnot(all(match(sub("_CHANGE-ME.*","",basename(dirname(filelist))),
#                    samples$SampleID) == 1:nrow(samples)))

#' * add filelist to samples as a new column
#samples %<>% mutate(Filenames = filelist)

#' * export full rank samples
#write_tsv(samples,here("doc/samples_full_rank.txt"))

#' Read the expression at the gene level
#' ```{r CHANGEME4,eval=FALSE,echo=FALSE}
#' If the species has only one transcript per gene, or if you are conducting QA 
#' in transcript level replace with the following:
#' txi <- suppressMessages(tximport(files = samples$Filenames, type = "salmon",txOut=TRUE))
#' ```
txi <- suppressMessages(tximport(files = samples$salmon_path,
                                 type = "salmon",
                                 tx2gene=tx2gene))
counts <- txi$counts
colnames(counts) <- samples$sample
#' 
#' <hr />
#' &nbsp;
#' 
#' # Quality Control
#' * "Not expressed" genes
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' ## Sequencing depth
#' * Let us take a look at the sequencing depth, colouring by CHANGEME
#' ```{r CHANGEME5,eval=FALSE,echo=FALSE}
#' # In the following most often you need to replace CHANGEME by your
#' # variable of interest, i.e. the metadata represented as column in
#' # your samples object, e.g. SamplingTime
#' ```
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

png(file.path(out_dir, "QA_no_merged_technical_reps", "sequencing_depth_sample_date.png"), width = 1200, height = 600)
ggplot(dat,aes(x,y,fill=Sampling_date)) + 
  geom_col() + 
  scale_y_continuous(name="reads") +
  facet_grid(~ factor(Sampling_date), scales = "free") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=4),
        axis.title.x=element_blank())
graphics.off()

png(file.path(out_dir, "QA_no_merged_technical_reps", "sequencing_depth_tissue.png"), width = 1200, height = 600)
ggplot(dat,aes(x,y,fill=Tissue_exp)) + 
  geom_col() + 
  scale_y_continuous(name="reads") +
  facet_grid(~ factor(Tissue_exp), scales = "free") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=4),
        axis.title.x=element_blank())
graphics.off()

#' `r emoji("point_right")` **We observe almost no difference in the raw sequencing depth**


#' ## per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
png(file.path(out_dir, "QA_no_merged_technical_reps", "per_gene_mean_expr.png"),   width = 1200, height = 600)
ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density(na.rm = TRUE) +
  ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)") + 
  theme_bw()
graphics.off()

#' `r emoji("point_right")` **The cumulative gene coverage is as expected**

#' ```{r CHANGEME6,eval=FALSE,echo=FALSE}
#' # In the following, the second mutate also needs changing, I kept it 
#' # as an example to illustrate the first line. SampleID would be 
#' # a column in the samples object (the metadata) that uniquely identify
#' # the samples.
#' # If you have only a single metadata, then remove the second mutate call
#' # If you have more, add them as needed.
#' ```
#' 
#' ## Per-sample expression
dat <- as.data.frame(log10(counts)) %>% 
  utils::stack() %>%
  mutate(Sampling_date=samples$Sampling_date[match(ind,samples$sample)]) %>% 
  mutate(Tissue_exp=samples$Tissue_exp[match(ind,samples$sample)])

png(file.path(out_dir, "QA_no_merged_technical_reps", "per_sample_expression_sample_date.png"), width = 1200, height = 600)
ggplot(dat,aes(x=values,group=ind,col=Sampling_date)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  theme_bw()
graphics.off()

png(file.path(out_dir, "QA_no_merged_technical_reps", "per_sample_expression_tissue.png"), width = 1200, height = 600)
ggplot(dat,aes(x=values,group=ind,col=Tissue_exp)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  theme_bw()
graphics.off()

#' `r emoji("point_right")` **All samples have the same sequencing depth characteristics and there is no deviation when we look at one or the other variable**
#' 
#' * Export raw expression data
write.csv(counts, file=file.path(out_dir, "QA_no_merged_technical_reps", "raw-unnormalised-gene-expression_data.csv"))
#' 
#' <hr />
#' &nbsp;
#' 
#' # Data normalisation 
#' ## Preparation
#' 
#' For visualization, the data is submitted to a variance stabilization
#' transformation using _DESeq2_. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
#'  ```{r CHANGEME7,eval=FALSE,echo=FALSE}
#'  # In the following, we provide the expected expression model, based on the study design.
#'  # It is technically irrelevant here, as we are only doing the quality assessment of the data, 
#'  # but it does not harm setting it correctly for the differential expression analyses that may follow.
#'  ```
dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData = samples,
  design = ~ 1) # no model 
colnames(dds) <- samples$sample

save(dds,file=file.path(out_dir, "QA_no_merged_technical_reps", "dds.rda"))

#' ## size factors 
#' (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds) %>% 
  suppressMessages()

png(file.path(out_dir, "QA_no_merged_technical_reps", "size_factors.png"),   width = 1200, height = 600)
boxplot(normalizationFactors(dds), # Took forever to load 
        main="Sequencing libraries size factor",
        las=2,log="y")
abline(h=1, col = "Red", lty = 3)
graphics.off()

#' and without outliers:
png(file.path(out_dir, "QA_no_merged_technical_reps", "size_factors_no_outlier.png"),  width = 1200, height = 600)
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor without outliers",
        las=2,log="y", outline=FALSE)
abline(h=1, col = "Red", lty = 3)
graphics.off()

#' Assess whether there might be a difference in library size linked to a
#' given metadata
#' ```{r echo=FALSE,eval=FALSE}
#' # Developer: This would need to be ggplot2'ed
#' ```
#boxplot(split(sizes,dds$Tissue),las=2,
#        main="Sequencing libraries size factor by Tissue")

#plot(sizes,log10(colSums(counts(dds))),ylab="log10 raw depth",xlab="scaling factor",
#     col=rainbow(n=nlevels(dds$Tissue))[as.integer(dds$Tissue)],pch=19)
#legend("bottomright",fill=rainbow(n=nlevels(dds$Tissue)),
#       legend=levels(dds$Tissue),cex=0.6)

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Validation
#' 
#' let's look at standard deviations before and after VST normalization. 
#' This plot is to see whether there is a dependency of SD on the mean. 
#' 
#' Before:  
png(file.path(out_dir, "QA_no_merged_technical_reps", "meanSDplot_dds.png"),  width = 1200, height = 800)
meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])
graphics.off()

#' After:
png(file.path(out_dir, "QA_no_merged_technical_reps", "meanSDplot_vst.png"),  width = 1200, height = 800)
meanSdPlot(vst[rowSums(vst)>0,])
graphics.off()

#' After VST normalization, the red line is almost horizontal which indicates
#' no dependency of variance on mean (homoscedastic).
#' 
#' `r emoji("point_right")` **We can conclude that the variance stabilization worked adequately**
#' 
#' <hr />
#' &nbsp;
#' 
#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' ### Scree plot
#' 
#' We define the number of variable of the model: ```r nvar=2```
#' `r nvar = 2`
nvar = 1 # just to set a nr 
#' An the number of possible combinations
#' ```{r CHANGEME8,eval=FALSE,echo=FALSE}
#' This needs to be adapted to your study design. Add or drop variables as needed.
#' ```
#' `r nlevel=nlevels(dds$MDay) * nlevels(dds$MGenotype)`
nlevel=1 # just to set a nr 
#' We plot the percentage explained by different components
#' 
#' * the red line represents number of variables in the model  
#' * the orange line represents number of variable combinations.
#' 
#' 
png(file.path(out_dir, "QA_no_merged_technical_reps", "Scree_plot.png"), width = 1200, height = 1200)
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)
graphics.off()

#' ### PCA plot
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    PC3=pc$x[,3],
                    as.data.frame(colData(dds)))

#PC1 vs PC2

pc.dat <- pc.dat %>% mutate(Genotype = as.factor(Genotype),
                            Geno_rep = paste0(Genotype, "_", Replicate))
dup_data <- pc.dat %>% select(c(Tissue_exp, Genotype, Replicate)) %>% filter(duplicated(.))
pc.dup_data <- left_join(dup_data, pc.dat, by =c("Tissue_exp", "Genotype", "Replicate"))

p1 <- ggplot(pc.dat,
             aes(x=PC1,y=PC2,col=Tissue_exp, shape=Genotype, text=Geno_rep)) + 
  geom_point(size=2) + 
  scale_shape_manual(values=1:nlevels(pc.dat$Genotype)) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") + 
  geom_point(data = pc.dup_data, aes(x = PC1, y = PC2, color = "red"), size = 5, shape = 21, show.legend = FALSE, alpha = 2/10)
  
my_plot <- ggplotly(p1) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep=""))) %>% suppressWarnings()
htmlwidgets::saveWidget(my_plot, file.path(out_dir, "QA_no_merged_technical_reps", "PCA_PC1_vs_PC2.html")) 

#PC1 vs PC3
p2 <- ggplot(pc.dat,aes(x=PC1,y=PC3,col=Tissue_exp, shape=Genotype, text=Geno_rep)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  geom_point(data = pc.dup_data, aes(x = PC1, y = PC3, color = "red"), size = 5, shape = 21, show.legend = FALSE, alpha = 2/10)

my_plot2 <- ggplotly(p2) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[3],"%)",sep=""))) %>% suppressWarnings()
htmlwidgets::saveWidget(my_plot2, file.path(out_dir, "QA_no_merged_technical_reps", "PCA_PC1_vs_PC3.html")) 

#' ```{r subplot, out.width = '100%'}
#' subplot(style(p1, showlegend = FALSE), p2,
#'         titleX = TRUE, titleY = TRUE, nrows = 1, margin = c(0.05, 0.05, 0, 0))
#' ```

#' ## Sequencing depth
#' Number of genes expressed per condition at different cutoffs:
conds <- factor(paste(dds$Tissue_exp,dds$Tissue_exp))
png(file.path(out_dir, "QA_no_merged_technical_reps", "Sequencing_depth.png"), width = 1200, height = 1200)
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)
graphics.off()

#' ## Heatmap
#' w
#' Filter for noise
#' 
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3) %>% suppressWarnings()
vst.cutoff <- 2

nn <- nrow(vst[sels[[vst.cutoff+1]],])
tn <- nrow(vst)
pn <- round(nn * 100/ tn, digits=1)

#' `r emoji("warning")` **`r pn`%** (`r nn`) of total `r tn` genes are plotted below:

#' ```{r CHANGEME9,eval=FALSE,echo=FALSE}
#' Optionally You can decide which variable to select for annotation. 
#' add the following to your pheatmap call (select the variable you want)
#' This will be a colored ribbon on top of the heatmap for better reading
#' annotation_col = samples %>% select(CHANGEME) %>% as.data.frame()
#' 
#' if you want, you can also modify the colors to better match with previous graphs. 
#' e.g. if you want a colored ribbon for Temperature (16,23,27) (adjust this example based on your need and add to pheatmap call)
#' annotation_colors = list(Temperature = c("23" = "#F8766D", "16" = "#00BF7D", "27" = "#00B0F6"))
#' 
#' Alternatively you can use the heatmap.2 version
#' #heatmap.2 version
#' hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
#'                 distfun=pearson.dist,
#'                 hclustfun=function(X){hclust(X,method="ward.D2")},
#'                 labRow = NA,
#'                 trace = "none",
#'                 labCol = conds,
#'                 col=hpal)
#' ```
#'  
mat <- t(scale(t(vst[sels[[vst.cutoff+1]],])))
png(file.path(out_dir, "QA_no_merged_technical_reps", "Heatmap.png"), width = 1200, height = 1200)
hm <- pheatmap(mat,
               color = hpal,
               border_color = NA,
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D2",
               show_rownames = FALSE,
               labels_col = conds,
               angle_col = 90,
               legend = FALSE)
graphics.off()

#' ## Clustering of samples
#' ```{r echo=FALSE,eval=FALSE}
#' # Developer: This wouldonly works with the gplots heatmap.2, not the pheatmap
#' plot(as.hclust(hm$colDendrogram),xlab="",sub="")
#' ```
#'
#' Below we assess the previous dendrogram's reproducibility and plot the clusters with au and bp where:
#' 
#' * __au (Approximately Unbiased): computed by multiscale bootstrap resampling__ `r emoji("point_left")` a better approximation
#' * __bp (Bootstrap Probability): computed by normal bootstrap resampling__
#' 
#' `r emoji("warning")`Clusters with AU larger than 95% are highlighted by rectangles, which are strongly supported by data
#' 
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                      method.hclust = "ward.D2", 
                      nboot = 50, parallel = TRUE) #nboot was around 1000 - I have too many samples to use that interactively 

png(file.path(out_dir, "QA_no_merged_technical_reps", "Clustering_of_samples.png"), width = 3000, height = 1200)
plot(hm.pvclust, labels = conds)
pvrect(hm.pvclust)
graphics.off()

#' <details><summary>bootstrapping results as a table</summary>
#' ```{r bootstrapping results as a table}
#' print(hm.pvclust, digits=3)
#' ```
#' </details>
#' 
#' ```{tech rep, echo=FALSE, eval=FALSE}
#' # The block of code is meant to combine tech reps - as it is facultative it is commented out
#' # First create a new variable in your sample object called BioID that identifies uniquely technical replicates, so one value for all tech rep of the same bio rep
#' samples$BioID <- CHANGEME
#' # Merging technical replicates
txi$counts <- sapply(split.data.frame(t(txi$counts),samples$Tube_nr),colSums)
txi$length <- sapply(split.data.frame(t(txi$length),samples$Tube_nr),colMaxs)
#' # Counts are now in alphabetic order, check and reorder if necessary
stopifnot(colnames(txi$counts) == samples$Tube_nr)
samples <- samples[match(colnames(txi$counts),samples$Tube_nr),]
#' # Recreate the dds
dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData = samples,
  design = ~ 1)
colnames(dds) <- samples$sample

# Save data
save(dds,file=file.path(out_dir, "dds.rda"))

dds <- estimateSizeFactors(dds) %>% 
  suppressMessages()

## size factors 
png(file.path(out_dir, "size_factors.png"),   width = 1200, height = 600)
boxplot(normalizationFactors(dds), # Took forever to load 
        main="Sequencing libraries size factor",
        las=2,log="y")
abline(h=1, col = "Red", lty = 3)
graphics.off()

#' and without outliers:
png(file.path(out_dir, "size_factors_no_outlier.png"),  width = 1200, height = 600)
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor without outliers",
        las=2,log="y", outline=FALSE)
abline(h=1, col = "Red", lty = 3)
graphics.off()

# vst 
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)
dim(vst) # 37075   215
vst <- vst[rowSums(vst)>0,]
dim(vst) # 33374   215

# Save data 
write.table(vst, file = file.path(out_dir, "exAtlas_vst.tsv"), sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)
write_tsv(colData(dds) %>% as.data.frame(), file = file.path(out_dir, "Meta_data_dds.tsv"))

## Validation 
#' Before:  
png(file.path(out_dir, "meanSDplot_dds.png"),  width = 1200, height = 800)
meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])
graphics.off()

#' After:
png(file.path(out_dir, "meanSDplot_vst.png"),  width = 1200, height = 800)
meanSdPlot(vst[rowSums(vst)>0,])
graphics.off()

# PCA 
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)
nvar = 1 
nlevel = 1
png(file.path(out_dir, "Scree_plot.png"), width = 1200, height = 1200)
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)
graphics.off()

#' ### PCA plot
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    PC3=pc$x[,3],
                    as.data.frame(colData(dds)))

#PC1 vs PC2
p1 <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Tissue_exp,shape=Tissue_exp,text=Tissue_exp)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

my_plot <- ggplotly(p1) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep=""))) %>% suppressWarnings()
htmlwidgets::saveWidget(my_plot, file.path(out_dir, "PCA_PC1_vs_PC2.html")) 

#PC1 vs PC3
p2 <- ggplot(pc.dat,aes(x=PC1,y=PC3,col=Tissue_exp,shape=Tissue_exp,text=Tissue_exp)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

my_plot2 <- ggplotly(p2) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[3],"%)",sep=""))) %>% suppressWarnings()
htmlwidgets::saveWidget(my_plot2, file.path(out_dir, "PCA_PC1_vs_PC3.html")) 

# Sequencing depth
conds <- factor(paste(dds$Tissue_exp,dds$Tissue_exp))
png(file.path(out_dir, "Sequencing_depth.png"), width = 1200, height = 1200)
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)
graphics.off()

## Heatmap
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3) %>% suppressWarnings()
vst.cutoff <- 2

nn <- nrow(vst[sels[[vst.cutoff+1]],])
tn <- nrow(vst)
pn <- round(nn * 100/ tn, digits=1)

mat <- t(scale(t(vst[sels[[vst.cutoff+1]],])))
png(file.path(out_dir, "Heatmap.png"), width = 2000, height = 1200)
hm <- pheatmap(mat,
               color = hpal,
               border_color = NA,
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D2",
               show_rownames = FALSE,
               labels_col = factor(paste(dds$Tissue_exp,dds$Genotype)),
               angle_col = 90,
               legend = FALSE)
graphics.off()

hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                      method.hclust = "ward.D2", 
                      nboot = 50, parallel = TRUE) #nboot was around 1000 - I have too many samples to use that interactively 

png(file.path(out_dir, "Clustering_of_samples.png"), width = 3200, height = 1200)
plot(hm.pvclust, labels = factor(paste(dds$Tissue_exp,dds$Genotype)))
pvrect(hm.pvclust)
graphics.off()

#'```
#'
#' <hr />
#' &nbsp;
#' 
#' # Summary
#' `r emoji("star")` **CHANGE-ME**
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' <details><summary>Session Info</summary>
#' ```{r session info}
#' sessionInfo()
#' ```
#' </details>
#'   
#' &nbsp;
#' 
