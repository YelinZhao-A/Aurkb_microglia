# QC and integration

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyseurat)

dir.home = getwd()
outdir = paste0(dir.home,"/output")
if (dir.exists(outdir) != 1) {
  dir.create(outdir)
}
setwd(outdir)

# pre-QC and QC each sample####
qc.out.dir <- paste0(outdir,"/QC_output")
if (dir.exists(qc.out.dir) != 1) {
  dir.create(qc.out.dir)
  print("create outdir")}

for (i in 1:length(all_sample.list)) {
  all_sample.list[[i]]$log10GenesPerUMI <- log10(all_sample.list[[i]]$nFeature_RNA)/log10(all_sample.list[[i]]$nCount_RNA)
  p <- VlnPlot(all_sample.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mito","log10GenesPerUMI"), ncol = 4)
  ggsave(plot=p,  path = qc.out.dir, dpi = 300, width = 45, height = 40, units = "cm" ,filename=paste0("QC-VlnPlot_",i,".pdf"))
  p <- FeatureScatter(object = all_sample.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(plot=p,  path = qc.out.dir, dpi = 300, width = 45, height = 40, units = "cm" ,filename=paste0("QC-feature_Plot1_",i,".pdf"))
  p <- FeatureScatter(object = all_sample.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mito")
  ggsave(plot=p,  path = qc.out.dir, dpi = 300, width = 45, height = 40, units = "cm" ,filename=paste0("QC-feature_Plot2_",i,".pdf"))
}

nCells_features <- matrix(NA,ncol = 3,nrow=length(all_sample.list))
for (i in 1:length(all_sample.list)) {
  nCells_features[i,1] = i
  nCells_features[i,2] = dim(all_sample.list[[i]])[1]
  nCells_features[i,3] = dim(all_sample.list[[i]])[2]
}
colnames(nCells_qc_features) = c("GSM","nGenes","nCells")
write.csv(nCells_features, paste0(outdir,"/nCells_features-preQC",".csv"),row.names = T, col.names = T, sep = ',',quote = F)

# QC each sample
all_sample_qc.list <- list()
for (i in 1:(length(all_sample.list))){ 
  # all_sample.list[i][['percent.mito']] <- PercentageFeatureSet(all_sample.list[i], pattern = "^MT-")
  all_sample.list[[i]]$log10GenesPerUMI <- log10(all_sample.list[[i]]$nFeature_RNA)/log10(all_sample.list[[i]]$nCount_RNA)
  all_sample_qc = subset(x = all_sample.list[[i]],
                         subset = percent.mito < 20 &
                           nFeature_RNA > 200 & nFeature_RNA < 6000 &
                           nCount_RNA > 400 &  nCount_RNA < 60000 &
                           log10GenesPerUMI > 0.80
  ) #Outlier cells in these samples might be cells that have a less complex RNA species than other cells. Sometimes we can detect contamination with low complexity cell types like red blood cells via this metric. Generally, we expect the novelty score to be above 0.80.
  all_sample_qc.list <- c(all_sample_qc.list,all_sample_qc)
}
head(all_sample_qc.list[[1]]@assays$RNA@counts)
all_sample_qc.list[[1]]@assays$RNA@counts[1:5,1:5]

nCells_qc_features <- matrix(NA,ncol = 3,nrow=length(all_sample_qc.list))
for (i in 1:length(all_sample_qc.list)) {
  nCells_qc_features[i,1] = i
  nCells_qc_features[i,2] = dim(all_sample_qc.list[[i]])[1]
  nCells_qc_features[i,3] = dim(all_sample_qc.list[[i]])[2]
}
colnames(nCells_qc_features) = c("GSM","nGenes","nCells")
write.csv(nCells_features, paste0(outdir,"/nCells_features-postQC",".csv"),row.names = T, col.names = T, sep = ',',quote = F)


save(all_sample_qc.list,file = paste0(outdir,"/all_sample_qc.list.RData"))

#Perform integration####
all_sample_qc.list <- lapply(X = all_sample_qc.list, FUN = function(x) {
  x <- NormalizeData(x) #Seurat的默认标准化方法是：每个细胞的某一count先除以该细胞的总count，然后乘以scale因子10000，再做个对数转换
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = all_sample_qc.list)
anchors <- FindIntegrationAnchors(object.list = all_sample_qc.list, 
                                  reference = c(1),
                                  anchor.features = features)
all.integrated <- IntegrateData(anchorset = anchors)
save(all.integrated,file = paste0(outdir,"/all_integrated.RData"))

#Analysis on integrated data####
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(all.integrated) <- "integrated"
seed = 42
nPCs = 15
res= 1
reduction= "tsne" 

standard10X = function(srat,nPCs=50,res=1.0,verbose=FALSE){
  srat = ScaleData(srat,verbose=verbose)
  srat = RunPCA(srat,verbose=verbose)
  p = ElbowPlot(srat,ndims = 30)
  ggsave(plot=p,  path = outdir, dpi = 300, width = 25, height = 20, units = "cm" ,filename="PCA-ElbowPlot.pdf")
  srat = RunTSNE(srat, seed.use = seed, reduction = "pca", dims=seq(nPCs))
  srat = FindNeighbors(srat, seed.use = seed, reduction = "pca", dims = seq(nPCs))
  srat = FindClusters(srat, random.seed = seed, resolution = res)
  return(srat)
}

all.integrated <- standard10X(all.integrated, nPCs=nPCs, res=res)

all.integrated@meta.data$ID = rownames(all.integrated@meta.data)

all.integrated@meta.data = merge(all.integrated@meta.data, meta.t, by.x = 'orig.ident', by.y="accession", all.x = T)
identical(all.integrated@meta.data$ID, colnames(all.integrated@assays$RNA))
rownames(all.integrated@meta.data) = all.integrated@meta.data$ID

save(all.integrated, file = paste0(outdir,"/all.integrated.RData"))
