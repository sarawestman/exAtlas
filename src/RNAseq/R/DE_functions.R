#  Plot specific gene expression
line_plot <- function(dds=dds,vst=vst,gene_id=gene_id){
  message(paste("Plotting",gene_id))
  sel <- grepl(gene_id,rownames(vst))
  stopifnot(sum(sel)==1)
  
  p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                        data.frame(value=vst[sel,])),
              aes(x=Treatment,y=value,col=Treatment,group=Treatment)) +
    geom_point() + geom_smooth() +
    scale_y_continuous(name="VST expression") + 
    ggtitle(label=paste("Expression for: ",gene_id))
  
  suppressMessages(suppressWarnings(plot(p)))
  return(NULL)
}

# Extract the DE results. Default cutoffs are from Schurch _et al._, RNA, 2016
extract_results <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("data/analysis/cold/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds),
                              expression_cutoff=0,
                              debug=FALSE,filter=c("median",NULL),...){
  
  # get the filter
  if(filter == "median"){
    filter2 <- rowMedians(counts(dds,normalized=TRUE))
    message("Using the median normalized counts as default, set filter=NULL to revert to using the mean")
  }
  
  # validation
  if(length(contrast)==1 & filter == "median"){
    res <- results(dds,name=contrast,filter = filter2)
  } else if(length(contrast)==1 & filter == "NULL"){
    res <- results(dds,name=contrast)
  } else {
    res <- results(dds,contrast=contrast,filter = filter2)
  }
  stopifnot(length(sample_sel)==ncol(vst))
  
  if(plot){
    pdf(file.path(default_dir, "Volcano_plot.pdf"))
    par(mar=c(5,5,5,5))
    volcanoPlot(res)
    par(mar=mar)
    graphics.off()
  }
  
  # a look at independent filtering
  if(plot){
    pdf(file.path(default_dir, "Independent_filtering_plot.pdf"))
    plot(metadata(res)$filterNumRej,
         type="b", ylab="number of rejections",
         xlab="quantiles of filter")
    lines(metadata(res)$lo.fit, col="red")
    abline(v=metadata(res)$filterTheta)
    graphics.off()
  }
  
  if(verbose){
    message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                    round(metadata(res)$filterThreshold,digits=5),
                    names(metadata(res)$filterThreshold)))
    
    max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
    message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                    round(max.theta*100,digits=5),
                    round(quantile(counts(dds,normalized=TRUE),probs=max.theta),digits=5)))
  }
  
  if(plot){
    qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
    dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                      qtl.exp=qtl.exp,
                      number.degs=sapply(lapply(qtl.exp,function(qe){
                        res$padj <= padj & abs(res$log2FoldChange) >= lfc & 
                          ! is.na(res$padj) & res$baseMean >= qe
                      }),sum))
    if(debug){
      plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
             geom_line() + geom_point() +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("base mean expression") +
             geom_hline(yintercept=expression_cutoff,
                        linetype="dotted",col="red"))
      
      p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
        geom_line() + geom_point() +
        scale_x_continuous("quantiles of expression") + 
        scale_y_log10("base mean expression") + 
        geom_hline(yintercept=expression_cutoff,
                   linetype="dotted",col="red")
      suppressMessages(suppressWarnings(plot(p)))
      
      plot(ggplot(dat,aes(x=thetas,y=number.degs)) + 
             geom_line() + geom_point() +
             geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("Number of DE genes"))
      
      plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) + 
             geom_line() + geom_point() +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("Cumulative number of DE genes"))
      
      plot(ggplot(data.frame(x=dat$thetas[-1],
                             y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
             geom_line() + geom_point() +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("Number of DE genes per interval"))
      
      plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                             y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
             geom_line() + geom_point() +
             scale_x_continuous("base mean of expression") + 
             scale_y_continuous("Number of DE genes per interval"))
      
      p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                             y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
        geom_line() + geom_point() +
        scale_x_log10("base mean of expression") + 
        scale_y_continuous("Number of DE genes per interval") + 
        geom_vline(xintercept=expression_cutoff,
                   linetype="dotted",col="red")
      suppressMessages(suppressWarnings(plot(p)))
    }
  }
  
  sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
    res$baseMean >= expression_cutoff
  
  if(verbose){
    message(sprintf(paste(
      ifelse(sum(sel)==1,
             "There is %s gene that is DE",
             "There are %s genes that are DE"),
      "with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s"),
      sum(sel),padj,
      lfc,expression_cutoff))
  }
  
  # proceed only if there are DE genes
  if(sum(sel) > 0){
    val <- rowSums(vst[sel,sample_sel,drop=FALSE])==0
    if (sum(val) >0){
      warning(sprintf(paste(
        ifelse(sum(val)==1,
               "There is %s DE gene that has",
               "There are %s DE genes that have"),
        "no vst expression in the selected samples"),sum(val)))
      sel[sel][val] <- FALSE
    } 
    
    if(export){
      if(!dir.exists(default_dir)){
        dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
      }
      write.csv(res,file=file.path(default_dir,paste0(default_prefix,"_results.csv")))
      write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"_genes.csv")))
    }
    pdf(file.path(default_dir, "Heatmap.pdf"), width = 30, height = 10)
    par(mar=c(0,5,5,5)+0.1)
    if(plot & sum(sel)>1){
      heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                distfun = pearson.dist,
                hclustfun = function(X){hclust(X,method="ward.D2")}, dendrogram='row',  Rowv=TRUE,
                trace="none",col=hpal,labRow = FALSE, srtCol = 45, margin=c(6,6),
                labCol=labels[sample_sel],...
      )
    }
    graphics.off()
  }
  return(list(all=rownames(res[sel,]),
              up=rownames(res[sel & res$log2FoldChange > 0,]),
              dn=rownames(res[sel & res$log2FoldChange < 0,])))
}

# Extract and plot the enrichment results
extractEnrichmentResults <- function(enrichment,task="go",
                                     diff.exp=c("all","up","dn"),
                                     go.namespace=c("BP","CC","MF"),
                                     genes=NULL,export=TRUE,plot=TRUE,
                                     default_dir=here("data/analysis/DE"),
                                     default_prefix="DE",
                                     url="athaliana"){
  # process args
  diff.exp <- match.arg(diff.exp)
  de <- ifelse(diff.exp=="all","none",
               ifelse(diff.exp=="dn","down",diff.exp))
  
  # sanity
  if( is.null(enrichment[[task]]) | length(enrichment[[task]]) == 0){
    message(paste("No enrichment for",task))
  } else {
    
    # write out
    if(export){
      write_tsv(enrichment[[task]],
                file=here(default_dir,
                          paste0(default_prefix,"-genes_GO-enrichment.tsv")))
      if(!is.null(genes)){
        write_tsv(
          enrichedTermToGenes(genes=genes,terms=enrichment[[task]]$id,url=url,mc.cores=16L),
          file=here(default_dir,
                    paste0(default_prefix,"-enriched-term-to-genes.tsv"))
        )
      }
    }
    
    if(plot){
      sapply(go.namespace,function(ns){
        titles <- c(BP="Biological Process",
                    CC="Cellular Component",
                    MF="Molecular Function")
        suppressWarnings(tryCatch({plotEnrichedTreemap(enrichment,enrichment=task,
                                                       namespace=ns,
                                                       de=de,title=paste(default_prefix,titles[ns]))},
                                  error = function(e) {
                                    message(paste("Treemap plot failed for",ns, 
                                                  "because of:",e))
                                    return(NULL)
                                  }))
      })
    }
  }
}