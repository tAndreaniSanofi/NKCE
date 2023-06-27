library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(pagoda2)
library(loomR)
library(tidyverse)
library(LoomExperiment)

main <- function(){
  #load velocyto loom file
  args <- commandArgs(trailingOnly = TRUE)
  loom_file <- args[1]
  marker_genes <- args[2]
  plot_figure <- args[3]
  gene_of_interest <- args[4]
  
  ldat <- ReadVelocity(loom_file)
  #extract differentially expressed genes in NKp46 or CTRL_Ab
  genes <- read.table(marker_genes)
  genesName <- genes$gene

  #Using spliced expression matrix as input to pagoda2.
  emat <- ldat$spliced
  length(rownames(emat))

  #this dataset has already been pre-filtered, but this is where one woudl do some filtering
  #emat <- emat[,colSums(emat)>=1e3]
  emat <- emat [rownames(emat) %in% genesName, ]
  length(rownames(emat))

  #Create pagoda2 object, adjust variance:
  rownames(emat) <- make.unique(rownames(emat))
  r <- Pagoda2$new(emat,modelType='plain',trim=20,log.scale=T)
  r$adjustVariance(plot=T,do.par=T,gam.k=10)

  #Run basic analysis steps to generate cell embedding and clustering, visualize:
  r$calculatePcaReduction(nPcs=100,n.odgenes=length(rownames(emat)),maxit=300)
  r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')
  r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
  r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)

  #Plot embedding, labeling clusters (left) and âFGFBP2â expression (which correspond to cluster 0 )

  par(mfrow=c(1,2))
  r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
  r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,"XCL1"],main='XCL1') 

  #Prepare matrices and clustering data
  emat <- ldat$spliced;
  dim(emat)
  emat <- emat [rownames(emat) %in% genesName, ]
  dim(emat)
  nmat <- ldat$unspliced;
  nmat <- nmat [rownames(nmat) %in% genesName, ]
  length(rownames(emat))
  length(rownames(nmat))

  #emat <- emat[,rownames(r$counts)]; # 
  #nmat <- nmat[,rownames(r$counts)]; # restrict to cells that passed p2 filter

  #take cluster labels
  cluster.label <- r$clusters$PCA[[1]]
  cell.colors <- sccore::fac2col(cluster.label)
  #take embedding
  emb <- r$embeddings$PCA$tSNE
  #in addition to clustering and the t-SNE embedding, from the p2 processing we will also take a cell-cell distance, which will be better than the default whole-transcriptome correlation distance that velocyto.R would normally use.
  cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))
  emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 1)
  nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.5)
  length(intersect(rownames(emat),rownames(nmat)))
  length(rownames(emat))
  length(rownames(nmat))
 
  #Estimate RNA velocity (using gene-relative model with k=20 cell kNN pooling and using top/bottom 2% quantiles for gamma fit):
  fit.quantile <- 0.02
  rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)

  #Visualize velocity on the t-SNE embedding, using velocity vector fields:
  show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)

    
}

intersect(rownames(emat),genesName)
gene <- "gene_of_interest"
gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 20,kGenes=1,fit.quantile=fit.quantile,cell.emb=emb,cell.colors=cell.colors,cell.dist=cell.dist,show.gene=gene,old.fit=rvel.cd,do.par=T)
intersect(rownames(emat),genesName)
print(Sample_8_annnotated_05_04)

#use the genes you prefer for cluster and trajectories annotation
bright_markers <- c( "FCGR3A", "GZMB", "CD56", "GZMK", "SELL", "CD2", "CCR7", "NCAM1")
CIML_markers <- c("XCL1", "XCL2", "KLRC1", "IL2RA")
typeI_inf <- c("HERC5", "IFIT2", "IFIT3", "OASL", "IFIT1", "ISG15", "ISG20", "PMAIP1", "ZC3HAV1", "RARRES3", "IFIH1", "IFI44", "RNF149", "PARP14", "CYTOR", "DRAP1", "RSAD2", "HERC6", "TMEM123", "SAR1A")
mature_markers <- c("KIR2DL1", "KIR3DL1", "KIR2DL3", "KIR3DL2", "S100A4", "ACTB", "ARPC3", "ARPC4", "CFL1", "PFN1", "CST7", "GNLY", "GZMB", "PRF1", "CD2")
