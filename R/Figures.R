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

##Figure 1a ####
load(file = paste0(dir.home,"/GSE121654/all.integrated.RData"))
p1 = DimPlot(all.integrated, reduction = "tsne", label = T)
p2 = DimPlot(all.integrated, reduction = "tsne", label = T, group.by = "age")
p = p1 + p2
ggsave(plot=p,  path = outdir, dpi = 300, width = 25, height = 12, units = "cm" ,filename="DimPlot_tsne.pdf")

#Feature plot
DefaultAssay(all.integrated) = 'RNA'
all.integrated$age = factor(all.integrated$age, levels=c('E14.5', 'P4','P5','P30','P100','P540' ))

Fig1a = FeaturePlot(all.integrated, features = 'Aurkb', max.cutoff = 3, label = FALSE,reduction = "tsne",
                     cols = c("grey", "red"), split.by = 'age') & theme(legend.position = "right")
ggsave(plot=Fig1a,  path = outdir, dpi = 300, width = 45, height = 6, units = "cm" ,filename=paste0("GSE121654_Feature_Plot_Aurkb_split.pdf"))


##Figure 1b ####
load(file = paste0(dir.home,"/GSE123025/all.integrated.RData"))
all.integrated$group = factor(all.integrated$group, levels=c('E14.5', 'P7','P60' ))
DefaultAssay(all.integrated) = 'RNA'
Fig1b = VlnPlot(
  all.integrated,
  features = c("Aurkb"),
  group.by = "group",
  pt.size = 0.1,  # Show individual cells
  log = TRUE      
) 
ggsave(filename = "/GSE123025_vlnplot_Aurkb.pdf", 
       plot = Fig1b, width = 3.5, height = 4)

##Figure 1d.1 ####
load(file = paste0(dir.home,"/GSE207570/all.integrated.RData"))
DefaultAssay(all.integrated) = 'RNA'
all.integrated$group = factor(all.integrated$group, levels=c('Naive', 'Demyelination'))
Fig1d = VlnPlot(
  all.integrated,
  features = c("Aurkb"),
  group.by = "group",
  pt.size = 0.1,  # Show individual cells
  log = TRUE      
)
ggsave(filename = "/GSE207570_vlnplot_Aurkb.pdf", 
       plot = Fig1b, width = 3.5, height = 4)

##Figure 1d.2 ####
load(file = paste0(dir.home,"/GSE204755/all.integrated.RData"))
DefaultAssay(all.integrated) = 'RNA'
all.integrated$group = factor(all.integrated$group, levels=c('Naive', 'Demyelination','Remyelination' ))
Fig1d2 = VlnPlot(
  all.integrated,
  features = c("Aurkb"),
  group.by = "group",
  pt.size = 0.1,  # Show individual cells
  log = TRUE      
)
ggsave(filename = "/GSE204755_vlnplot_Aurkb.pdf", 
       plot = Fig1d2, width = 4, height = 4)

##Figure 1f ####
load(file = paste0(dir.home,"/GSE301908/all.integrated.RData"))
DefaultAssay(all.integrated) = 'RNA'
all.integrated$group = factor(all.integrated$group, levels=c('Control', 'MS' ))
Fig1f = VlnPlot(
  all.integrated,
  features = c("AURKB","MKI67"),
  group.by = "group",
  pt.size = 0.1,  # Show individual cells
  log = TRUE      
)
ggsave(filename = "/GSE301908_vlnplot_Aurkb_Mki67.pdf", 
       plot = Fig1f, width = 7, height = 4)

##Figure S1a ####
load(file = paste0(dir.home,"/GSE123025/all.integrated.RData"))
DefaultAssay(all.integrated) = 'RNA'
all.integrated$group = factor(all.integrated$group, levels=c('E14.5', 'P7','P60' ))
FigS1a_aurkb = FeaturePlot(all.integrated, features = 'Aurkb', max.cutoff = 3, label = FALSE,reduction = "umap",
                     cols = c("grey", "red"), split.by = 'group') & theme(legend.position = "right")
FigS1a_mki67 = FeaturePlot(all.integrated, features = 'Mki67', max.cutoff = 3, label = FALSE,reduction = "umap",
                     cols = c("grey", "red"), split.by = 'group') & theme(legend.position = "right")
FigS1a = FigS1a_aurkb + FigS1a_mki67
ggsave(plot=FigS1a,  path = outdir, dpi = 300, width = 45, height = 6, units = "cm" ,filename=paste0("GSE123025_Feature_Plot_Aurkb_Mki67_split.pdf"))

##Figure S1b ####
load(file = paste0(dir.home,"/GSE301908/all.integrated.RData"))
DefaultAssay(all.integrated) = 'RNA'
all.integrated$group = factor(all.integrated$group, levels=c('Control', 'MS' ))
p = FeatureScatter(all.integrated, feature1 = "AURKB", feature2 = "MKI67", group.by = "group")
ggsave(filename = paste0(outdir, "/GSE301908_FeatureScatter_AURKB_MKI67.pdf"), 
       plot = p, width = 6, height = 4)


