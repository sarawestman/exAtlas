# https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html#Module_Construction
suppressPackageStartupMessages({
  library(WGCNA)
  library(flashClust)
  library(curl)
  library(tidyverse)
  library(RColorBrewer)
  source("~/src/UPSCb-common/src/R/gopher.R")
  source("~/Git/exAtlas/src/RNAseq/R/plotEnrichedTreemap.R")
})

# Set colors
my_palette = c(brewer.pal(11, "RdBu")[c(10,6,2)])
myColor <- colorRampPalette(my_palette)(100)
myBreaks <- c(seq(-1, 0, length.out = ceiling(100/2) + 1), # Use floor and ceiling to deal with even/odd length pallette lengths
              seq(1/100, 1, length.out = floor(100/2)))
# Load data 
sel_tissues <- c("Leaf_D3_GH")

out_dir <- file.path("~/Sara/exAtlas/res/RNAseq/analysis/WGCNA", paste(sel_tissues, collapse = "_"))
dir.create(out_dir, showWarnings = FALSE)

SPG_known_DAM <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/RC_quantification/SPG_DAM_singletons_rm.csv", sep = ";")
exAtlas_LCMS_RNA_meta_info <- read.delim("~/Sara/exAtlas/res/RNAseq/analysis/DE/exAtlas_LCMS_RNA_meta_info.tsv") %>% 
  filter(Tissue_exp %in% sel_tissues)
LCMS_dat <- read.delim("~/Sara/exAtlas/res/Preprocessing_res/MS/DEM/Tissue_exp/exAtlas_LCMS_vsn_normalized_intensity.tsv")
LCMS_untarg_sel <- left_join(exAtlas_LCMS_RNA_meta_info, LCMS_dat, by = c(LCMS_sample = "Sample"))
expression.data <- read.delim("~/Sara/exAtlas/res/RNAseq/analysis/salmon/exAtlas_vst.tsv") %>% 
  `colnames<-`(gsub("\\_.*","", colnames(.))) %>% 
  as.data.frame() %>% 
  dplyr::select(all_of(exAtlas_LCMS_RNA_meta_info %>% dplyr::filter(Tissue_exp %in% sel_tissues) %>% pull(RNA_sample))) %>% 
  t() #transforming the data.frame so columns now represent genes and rows represent samples

#Identifying Outlier Genes
#The WGCNA package has a built in function to identify outlier genes called goodSampleGenes(). The function checks the data and returns a list object of samples and genes that pass its filtering criteria. You can adjust how strict the filtering process is by changing the default of the arguments (which can be found under the goodSampleGene documentation, ?goodSampleGenes).
gsg <- goodSamplesGenes(expression.data)

#By viewing the gsg list object you can see it contains 3 logical vectors (good genes, good samples and allOK). If you want to see if the function identified any possible outlier all you have to do is evaluate the allOK vector.
summary(gsg)
gsg$allOK
# If the allOK object returns true, which it does in this case, there are no outliers present. If it doesn’t return true, we need to filter out the outliers manually using the following code.

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(expression.data)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
  expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}

# Identifying outlier samples
# You can identify outlier samples by using hierarchical clustering. It is probably beneficial to review the clustering module as clustering is used several times within this pathway.
sampleTree <- hclust(dist(expression.data), method = "average") #Clustering samples based on distance 
#Setting the graphical parameters
par(cex = 0.6)
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# You can remove the outlier using a cutree function.
#Setting the graphical parameters
#pdf(file.path(out_dir, "Sample_outlier.pdf"))
#par(cex = 0.6)
#par(mar = c(0,4,2,0))
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#     cex.axis = 1.5, cex.main = 2)
#draw on line to show cutoff height
#abline(h = 130, col = "red")
#graphics.off()

#cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 130, minSize = 10) #returns numeric vector
#Remove outlier
#expression.data <- expression.data[cut.sampleTree==1, ]

# Network Construction
spt <- pickSoftThreshold(expression.data) 

# You can then plot this data frame to better visualize what β value you should choose.
# Plot the R2 values as a function of the soft thresholds
pdf(file.path(out_dir, "Soft_threshold_R2.pdf"))
par(mar = c(1, 1, 1, 1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2], col = "red")
abline(h = 0.80, col = "red")
graphics.off()

# Plot mean connectivity as a function of soft thresholds
pdf(file.path(out_dir, "Soft_threshold_mean_connectivity.pdf"))
par(mar = c(1, 1, 1, 1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels = spt$fitIndices[,1], col = "red")
graphics.off()
# You can determine the soft power threshold should be set to 6 as it is the spt that retains the highest mean connectivity while reaching an R2 value above 0.80.

# Calling the Adjacency Function
# Now that you have the soft threshold power determined you can call on the adjacency() function of the WGCNA package.
softPower <- 6
adjacency <- adjacency(expression.data, power = softPower)

# Module Construction
# Topological Overlap Matrix
# To convert the adjacency matrix into a TOM similarity matrix we can call the WGCNA function TOMsimilarity().
TOM <- TOMsimilarity(adjacency)

# To convert this matrix into a dissimilarity matrix you can subtract the TOM object from 1.
TOM.dissimilarity <- 1-TOM

# Hierarchical Clustering Analysis
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 

#plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)

# To identify modules from this gene dendrogram, you can use the cutreeDynamic() function. This will allow you to set a minimum cluster size. 
# For genomic data like this it is more beneficial to set minimum module sizes relatively high as you are working with high loads of data. 
# The authors of WGCNA recommend to start at a minClusterSize = 30.
Modules <- cutreeDynamic(dendro = geneTree, 
                         distM = TOM.dissimilarity, 
                         deepSplit = 2, 
                         pamRespectsDendro = FALSE, 
                         minClusterSize = 30)
table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module. 

#You can now plot the module assignment under the gene dendrogram for visualization.
ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors)

# find Hub genes 
#res <- chooseTopHubInEachModule(expression.data, ModuleColors) # in each module - looking at all genes in GE data 
#res <- chooseOneHubInEachModule(expression.data, ModuleColors) # in each module - given a random set of genes 

#plots the gene dendrogram with the module colors
pdf(file.path(out_dir,"TOM_based_dissimilarity_colors.pdf"))
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
graphics.off()

# Module Eigengene Identification
# A ME (Module Eigengene) is the standardized gene expression profile for a given module.
# To identify the Module Eigengene you can call on the expression data into the moduleEigengenes() function.
MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)
# You have now identified the eigengenes for the data.

# Module Merging
# To further condense the clusters (branches) into more meaningful modules you can cluster modules based on pairwise eigengene correlations and merge the modules that have similar expression profiles.
# In order to do cluster analysis you must first find a measurement of dissimilarity (distance) between module eigengenes.
# However, because there are missing values present in the input variable, you will need to add an argument use= "complete". This command removes rows of the matrix which have NA values. Removing these NAs allows ME.dissimilarity to run. This may or may not be necessary depending on your data set.
#my_dat <- MElist$eigengenes %>% as.data.frame() %>% mutate_all(~ifelse(is.nan(.), NA, .)) %>% select(-c(MEgrey)) # for Leaf_D1
#ME.dissimilarity = 1-cor(my_dat, use="complete") #Calculate eigengene dissimilarity
ME.dissimilarity = 1-cor(MElist$eigengenes, use = "complete") #Calculate eigengene dissimilarity

# Now, using the new found measurements of dissimilarity, you can construct a cluster tree. You will also be adding a line at the height of .25. This height corresponds to a correlation of over 75%. Any branches below this line are more than 75% related, and you will thus be merging them
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(METree)
abline(h = .25, col = "red") #a height of .25 corresponds to correlation of .75
#This figure shows all of the modules which are more than 75% similar. For example you can see that MEcyan and MEpurple are more than 75% similar. Now you can merge these two modules, and others like them using the mergeCloseModules() command

merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)
# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

pdf(file.path(out_dir,"TOM_based_dissimilarity_merged_modules.pdf"))
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
graphics.off()
# Coooollioo! you can see that generally, there are less colors in the merged modules row, showing that our modules which were more than 75% similar have merged. For example, the lime green and red section merged into one large lime green section. This new merged output is cleaner and further identifies groups of highly correlated genes which co-occur across samples.

# Save results 
save(geneTree, TOM.dissimilarity, expression.data, Modules, MElist, adjacency, ModuleColors, merge, file = file.path(out_dir, paste0("exAtlas_WGCNA_softpower",softPower,".RData")))
# ModuleColors will be missing for individual tissues (ie. Shoot_tip_GH, etc.)

Module_res <- data.frame(Genes = colnames(expression.data), modules = mergedColors)
write_tsv(Module_res, file = file.path(out_dir, "TOM_modules.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# External Trait Matching
# Once you have constructed the network (the adjacency matrix) and divided it into meaningful modules, you can begin to relate the network to external traits.

# You can clean the data by removing unnecessary columns and pulling out the continuous traits.
allTraits <- LCMS_untarg_sel %>% 
  column_to_rownames("RNA_sample") %>% 
  select(-c(LCMS_sample:isRoot)) %>% 
  `colnames<-`(str_extract(colnames(.), "[^_]+")) %>% 
  dplyr::select(all_of(paste0("RC", SPG_known_DAM$RC))) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("RC") %>% 
  left_join(., SPG_known_DAM %>% dplyr::select(-c(SA:Note)) %>% mutate(RC = paste0("RC", RC)), by = "RC") %>% 
  mutate(lab = paste0(Trait, " (", RC, ")")) %>% 
  dplyr::select(-c(RC, Trait)) %>%
  column_to_rownames("lab") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample")

my_cor <- cor(allTraits %>% dplyr::select(-c(Sample)))

# Add annotation 
annoRow <- SPG_known_DAM %>% 
  mutate(lab = ifelse(!is.na(Trait), paste0(Trait," (RC", RC, ")"), RC)) %>% 
  arrange(match(lab, colnames(allTraits))) %>% 
  column_to_rownames("lab") %>% 
  dplyr::select(-c(RC, Trait, Anno, Note, MZ, RT))
stopifnot(rownames(annoRow) == colnames(allTraits[,-1]))

ann_colors = list(
  CN_DAM = c(Up = "red", Dn = "blue", NS = "gray"),
  SA = c(Yes = "coral3", No = "gray"), 
  HCH = c(Yes = "coral3", No = "gray"),
  CN = c(Yes = "coral3", No = "gray"),
  BA = c(Yes = "coral3", No = "gray"),
  AC = c(Yes = "coral3", No = "gray"))

ph1 <- pheatmap::pheatmap(my_cor,
                          color = myColor,
                          annotation_row = annoRow,
                          annotation_colors = ann_colors,
                          breaks = myBreaks)

pdf(file.path(out_dir, "SPG_pot_heatmap.pdf"), width = 17, height = 15)
ph1
graphics.off()

allTraits <- allTraits[,c("Sample", rownames(my_cor[ph1$tree_row[["order"]],]))]
#write_tsv(allTraits, file.path(out_dir, "SPGs_allTraits.tsv"))

# Finally, you must match the trait data to the expression data by the sample number
Samples <- rownames(expression.data)
traitRows <- match(Samples, allTraits$Sample)
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]

# Module-Trait associations
# First, you will quantify the association between the expression profile and a particular trait of interest by calculating the correlation of the trait with previously identified module eigengenes. This pairwise correlation is known as the eigengene gene significance

# Define numbers of genes and samples
nGenes = ncol(expression.data)
nSamples = nrow(expression.data)
mergedMEs = merge$newMEs
module.trait.correlation = cor(mergedMEs, datTraits, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation
write_tsv(module.trait.correlation %>% as.data.frame() %>% rownames_to_column("modules"), 
          file.path(out_dir, "Module_SPGtrait_correlation.tsv"))

# Once you have the gene signficance and the corresponding p-value for all the modules and traits, you can create a graphical representation (module-trait heatmap) that will be helpful to neatly visualize the module-trait relationships.
# Will display correlations and their p-values
#module.trait.correlation <- module.trait.correlation[-53,] # remove MEgrey from Leaf_D1
#module.trait.Pvalue <- module.trait.Pvalue[-53,] # remove MEgrey from Leaf_D1
#mergedMEs <- mergedMEs[,-53] # remove MEgrey from Leaf_D1
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
                   signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
# Display the correlation values within a heatmap plot

pdf(file.path(out_dir,"TOM_SPG_pot.pdf"), height = 20, width = 20)
par(mar = c(12, 12, 10, 12))
labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
graphics.off()

pdf(file.path(out_dir,"TOM_SPG_pot_pheatmap.pdf"), height = 20, width = 20)
par(mar = c(12, 12, 10, 12))
pheatmap::pheatmap(module.trait.correlation,
                   color = myColor,
                   annotation_col = annoRow,
                   annotation_colors = ann_colors,
                   breaks = myBreaks,
                   fontsize = 8,
                   display_numbers = textMatrix)
graphics.off()

pdf(file.path(out_dir,"TOM_SPG_pot_pheatmap_noText.pdf"), height = 20, width = 20)
par(mar = c(12, 12, 10, 12))
pheatmap::pheatmap(module.trait.correlation,
                   color = myColor,
                   annotation_col = annoRow,
                   annotation_colors = ann_colors,
                   breaks = myBreaks)
graphics.off()

# Each row corresponds to a module eigengene, and the columns correspond to a trait. Each cell contains a p-value and correlation. Those with strong positive correlations are shaded a darker red while those with stronger negative correlations become more blue.

# Target Gene Identification
# Define variable weight containing the weight column of datTrait
weight_RC <- "HCH-cinnamoyl-salicortin_peak2 (RC465)"
weight = as.data.frame(datTraits[[weight_RC]])
names(weight) = "weight"
modNames = substring(names(mergedMEs), 3) #extract module names

#Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(cor(expression.data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")
write_tsv(MMPvalue, file = file.path(out_dir, "Module_membership_pval.tsv"))

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(cor(expression.data, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(weight), sep = "")
names(GSPvalue) = paste("p.GS.", names(weight), sep = "")
head(GSPvalue)
write_tsv(GSPvalue %>% arrange(p.GS.weight) %>% rownames_to_column("Gene"), file = file.path(out_dir, paste0("Gene_significance_pval_",str_replace_all(weight_RC, fixed(" "), ""),".tsv")))
# Using the gene significance you can identify genes that have a high significance for weight. Using the module membership measures you can identify genes with high module membership in interesting modules.
