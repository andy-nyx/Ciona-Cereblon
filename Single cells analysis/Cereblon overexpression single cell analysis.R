library("hdf5r")
library(cowplot)
library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(limma)
library(data.table)
library(stringr)
library(plyr)
library(tidyverse)
library(tidyft)
library("ComplexHeatmap")
#library("clustifier")

library("limma")
library(data.table)
library("RColorBrewer")
library(scico)

# Plotting
library(ggraph)
library("clustree")
library(viridisLite)
library("viridis")
library(circlize)

setwd("/your_directory") # Set Working directory
object <- "yourRobject.RData" #File name of the R object
ids <- c("Tbx6_GFP", "Tbx6_CRBN")

wdGFP <-"/controlDirectory" #Directory of the gene expression matrix for the control condition
wdCrbn <- "/cerblonOverexpressionDirectory" #Directory of the gene expression matrix for the Cereblon overexpression condition
filterMatrix <- "/filtered_feature_bc_matrix.h5"

GFP <- Read10X_h5(paste(wdGFP,filterMatrix, sep = ""), use.names = TRUE, unique.features = TRUE)
colnames(GFP) <- paste(sapply(strsplit(colnames(GFP),split="-"),'[[',1L),"Tbx6_GFP",sep="-")

Crbn <- Read10X_h5(paste(wdCrbn,filterMatrix, sep = ""), use.names = TRUE, unique.features = TRUE)
colnames(Crbn) <- paste(sapply(strsplit(colnames(Crbn),split="-"),'[[',1L),"Tbx6_CRBN",sep="-")

d10x.data <- list(GFP, Crbn)

names(d10x.data) <- ids

experiment.data <- do.call("cbind", d10x.data)


gfpvscrbn <- CreateSeuratObject(counts =  experiment.data, 
                                project = "GFPvsCRBN", 
                                min.cells = 3, 
                                min.features = 200,
                                names.field = 2,
                                names.delim = "\\-")

save.image(object)

pdf(file="gene-RNA-raw.content-1.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
VlnPlot(gfpvscrbn, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

pdf(file="gene-RNA-raw.content-2.pdf", 
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeatureScatter(gfpvscrbn, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()


# subset: cells with more than 2000 genes (nfeature_RNA)
# and outliners with high mRNA content (>70000 nCount_RNA)
gfpvscrbn <- subset(gfpvscrbn, subset = nFeature_RNA > 2000 
            & nCount_RNA < 70000)

pdf(file="gene-RNA.content-1.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
VlnPlot(gfpvscrbn, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) 
dev.off()

pdf(file="gene-RNA.content-2.pdf", 
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeatureScatter(gfpvscrbn, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# split the dataset into a list of two seurat objects

gfpvscrbn.list <- SplitObject(gfpvscrbn, split.by = "orig.ident")


# normalize and identify variable features for each dataset independently

gfpvscrbn.list <- lapply(X = gfpvscrbn.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration

features <- SelectIntegrationFeatures(object.list = gfpvscrbn.list)

#Perform integration

gfpvscrbn.anchors <- FindIntegrationAnchors(object.list = gfpvscrbn.list, anchor.features = features)

# this command creates an 'integrated' data assay

gfpvscrbn.combined <- IntegrateData(anchorset = gfpvscrbn.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay

DefaultAssay(gfpvscrbn.combined) <- "integrated"


# Run the standard workflow for visualization and clustering

gfpvscrbn.combined <- ScaleData(gfpvscrbn.combined, verbose = FALSE)
gfpvscrbn.combined <- RunPCA(gfpvscrbn.combined, npcs = 40, verbose = FALSE)

print(gfpvscrbn.combined[["pca"]], dims = 1:5, nfeatures = 5)
#DimHeatmap(gfpvscrbn.combined, dims = 1:10, cells = 500, balanced = TRUE)

pdf(file="PCA.pdf", 
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(gfpvscrbn.combined, reduction = "pca")
dev.off()

#Choosing number of PCAs for clustering

pdf(file="Elbowplot.pdf", 
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
ElbowPlot(gfpvscrbn.combined, ndims = 40)
dev.off()

#number of PCA for the analysis 20  (From 23.03.25 analysis_CCA, the resolution was too high, decrease the number of PCAs from 24 to 20)

n.dims <- 20

gfpvscrbn.combined <- FindNeighbors(gfpvscrbn.combined, reduction = "pca", dims = 1:n.dims)
resolutions <- seq(0.5, 1, 0.1)
gfpvscrbn.combined <- FindClusters(gfpvscrbn.combined,
                                   reduction.type = "pca",
                                   dims.use = 1:n.dims, 
                                   resolution = resolutions)


pdf(file="cluster-tree.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
clustree(gfpvscrbn.combined)
dev.off()


save.image(object)

# take resolution 0.7
head(Idents(gfpvscrbn.combined))

gfpvscrbn.combined$seurat_clusters <- gfpvscrbn.combined$integrated_snn_res.0.7
Idents(gfpvscrbn.combined) <- "seurat_clusters"

gfpvscrbn.combined <- RunTSNE(gfpvscrbn.combined, reduction = "pca", dims = 1:n.dims)

gfpvscrbn.combined <- RunUMAP(gfpvscrbn.combined, reduction = "pca",
                              dims = 1:n.dims, reduction.name = "UMAP")

gfpvscrbn.combined$orig.ident <- factor(x= gfpvscrbn.combined$orig.ident, levels = c("Tbx6_GFP", "Tbx6_CRBN"))

col1 <- c("#35608D","#66CC33") #Blue for control and Green for Crbn

save.image(file = object)

#visualization using t-SNE

pdf(file="sampleID-tSNE_.pdf", 
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(gfpvscrbn.combined, reduction = "tsne", 
        group.by = "orig.ident", cols = alpha(col1, 0.5))+
  coord_fixed()
dev.off()

DimPlot(gfpvscrbn.combined, reduction = "tsne", 
        group.by = "orig.ident", cols = alpha(col1, 0.5), label.size = 2, pt.size = 0.2)+
  coord_fixed()
ggsave("sampleID-tSNE-1d.pdf", device= "pdf", width = 20, 
       height = 16, units = "cm") 


DimPlot(gfpvscrbn.combined, reduction = "tsne", label = TRUE) + NoLegend()
ggsave("tSNE-cluster_res0.7.pdf", device= "pdf", width = 20, 
       height = 16, units = "cm") 

pdf(file="UMAP_seurat-cluster_res0.7.pdf", 
    width = 9, # The width of the plot in inches
    height = 8) # The height of the plot in inches
DimPlot(gfpvscrbn.combined, reduction = "UMAP", label = TRUE) + NoLegend()
dev.off()

DimPlot(gfpvscrbn.combined, reduction = "tsne", group.by = "integrated_snn_res.0.7", label = TRUE) + NoLegend()
DimPlot(gfpvscrbn.combined, reduction = "UMAP", group.by = "integrated_snn_res.0.6", label = TRUE) + NoLegend()

save.image(object)


# To visualize gene expression it should be done on the RNA assay and not integrated

col2 = c("lightgrey", "#330066")

DefaultAssay(gfpvscrbn.combined) <- "RNA"

pdf(file="KH2012:KH.C7.285-tSNE.pdf", #Cholinergic receptor
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(gfpvscrbn.combined, features = "KH2012:KH.C7.285",
            split.by = "orig.ident", reduction = "tsne",
            min.cutoff = "q10", max.cutoff = "q90", cols = col2)
dev.off()

pdf(file="KH2012:KH.C1.1116-tSNE.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(gfpvscrbn.combined, features = "KH2012:KH.C1.1116", #Hand2
            split.by = "orig.ident", cols = col2, reduction = "tsne",
            min.cutoff = "q10", max.cutoff = "q90")
dev.off()

pdf(file="KH2012:KH.C14.307-tSNE.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(gfpvscrbn.combined, features = "KH2012:KH.C14.307", #Myod
            reduction = "tsne",
            split.by = "orig.ident", cols = col2)
dev.off()

FeaturePlot(gfpvscrbn.combined, features = "KH2012:KH.C14.307",
            min.cutoff = "q10", max.cutoff = "q90", reduction = "tsne",
            label.size = 2, pt.size = 0.2, cols = col2) + coord_fixed()
ggsave("KH2012:KH.C14.307_MyoD-tSNE_2.pdf", device= "pdf", width = 20, 
       height = 16, units = "cm") 


pdf(file="CFPSV40-tSNE.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(gfpvscrbn.combined, features = "CFPSV40",
            split.by = "orig.ident", cols = col2, reduction = "tsne",
            min.cutoff = "q10", max.cutoff = "q90")
dev.off()

save.image(file = object)

# Identify cluster types (ATTENTION: as to be performed on the RNA assay slot)

DefaultAssay(gfpvscrbn.combined) <- "RNA"  
gfpvscrbn.combined.markers <- FindAllMarkers(gfpvscrbn.combined, only.pos = TRUE,
                                             min.pct = 0.4,
                                             logfc.threshold = 0.25,
                                             test.use = "roc",
                                             slot = "data")

#annotate markers
human.homo <- read.table("geneModel.txt", header = TRUE, row.names = 1)
human.homo <- human.homo[,c(1,2)]
colnames(human.homo)<-c("KH.model.ID","Human.homolog")
gfpvscrbn.combined.markers <- tidyft::left_join(gfpvscrbn.combined.markers,human.homo,
                                                by = c("gene" = "KH.model.ID"))
write.table(gfpvscrbn.combined.markers, file="DEG_44clusters_combined_res0.7.txt",
            quote=F, sep="\t", col.names=NA)

save.image(file = object)


#Germ cells
pdf(file="Germ cells-tSNE.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(gfpvscrbn.combined, features = c("KH2012:KH.C11.383",
                                             "KH2012:KH.C1.755",
                                             "KH2012:KH.C8.370",
                                             "KH2012:KH.S852.2"),
            reduction = "tsne", cols = col2,
            min.cutoff = "q10", max.cutoff = "q90")
dev.off()

#Prop cells

pdf(file="Eminens neurons-tSNE.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(gfpvscrbn.combined, features = c("KH2012:KH.C11.543", "KH2012:KH.C3.474"),
            reduction = "tsne", cols = col2,
            min.cutoff = "q10", max.cutoff = "q90")
dev.off()

#Pigment cells
pdf(file="Pigment cells-tSNE.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(gfpvscrbn.combined, features = c("KH2012:KH.C12.469",
                                             "KH2012:KH.C5.485",
                                             "KH2012:KH.C8.537",
                                             "KH2012:KH.C10.106"),
            reduction = "tsne", cols = col2,
            min.cutoff = "q10", max.cutoff = "q90")
dev.off()

pdf(file="BTNs-tSNE.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(gfpvscrbn.combined, features = c("KH2012:KH.C2.42", "KH2012:KH.C1.215"),
            reduction = "tsne", cols = col2,
            min.cutoff = "q10", max.cutoff = "q90")
dev.off()

pdf(file="Switch neurons-tSNE.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(gfpvscrbn.combined, features = c("KH2012:KH.C2.560", "KH2012:KH.L147.32"),
            reduction = "tsne", cols = col2,
            min.cutoff = "q10", max.cutoff = "q90")
dev.off()

pdf(file="VP neurons-tSNE.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(gfpvscrbn.combined, features = c("KH2012:KH.C6.11",
                                             "KH2012:KH.C3.635", 
                                             "KH2012:KH.C11.631"),
            reduction = "tsne", cols = col2,
            min.cutoff = "q10", max.cutoff = "q90")
dev.off()


save.image(file = object)

# Update annotation for cell types and tissues (use same annotation as 23.03.26 CCA analysis)

tissue.type2 <- read.csv("tissue.annotation.csv",header=TRUE)
tissue.type2 <- subset(tissue.type2, select = c("cluster","cell.type.3","tissue.2"))
colnames(tissue.type2)[2] <-"cell.type.2"

head(gfpvscrbn.combined@meta.data)

metadata <- gfpvscrbn.combined@meta.data
head(metadata)

md3 <- subset(metadata, select = seurat_clusters)

md3 <- rownames_to_column(md3)

head(md3)

class(md3$seurat_clusters)
class(tissue.type2$cluster)
md3[,"seurat_clusters"]<-as.character(md3[,"seurat_clusters"])
md3[,"seurat_clusters"]<-as.integer(md3[,"seurat_clusters"])

md3 <- left_join(md3,tissue.type2, by = c("seurat_clusters" = "cluster"))

head(md3)

md3 <- column_to_rownames(md3, var = "rowname")

head(md3)

md3 <- subset(md3,select = c(tissue.2, cell.type.2))

gfpvscrbn.combined <- AddMetaData(object = gfpvscrbn.combined,
                                  metadata = md3)


head(gfpvscrbn.combined@meta.data)

gfpvscrbn.combined$tissue.2<- factor(x= gfpvscrbn.combined$tissue.2, levels = c("epidermis",
                                                                             "nervous system",
                                                                             "endoderm",
                                                                             "mesenchyme",
                                                                             "trunk ventral cells",
                                                                             "muscle",
                                                                             "notochord",
                                                                             "germ cells"))


#Paired palette

nb.cols <- 44
col4a <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
barplot(1:44, col = col4a)

col5epi <- col4a[2]       #Blue
col5NS <- col4a[11]       #Green
col5endo <- "#EB4445" #Red
col5mes <- col4a[26] #Orange
col5TVCmusc <- c("#CAB2D6", "#6A3D9A") #Violet
col5germ<- col4a[43] #Brown
col5noto <- col4a[41] #Yellow
      
col5 <- c(col5epi, col5NS, col5endo, col5mes, col5TVCmusc, col5noto, col5germ)

barplot(1:8, col = col5)
        
DimPlot(gfpvscrbn.combined, reduction = "UMAP", 
                group.by = "tissue.2",
                label = FALSE,  label.size = 2, pt.size = 0.2, cols = col5)+
          coord_fixed()
ggsave("tissue.2-UMAP_paired.pdf", device= "pdf", width = 20, 
               height = 16, units = "cm") 
        
DimPlot(gfpvscrbn.combined, reduction = "tsne", 
                group.by = "tissue.2",
                label = FALSE,  label.size = 2, pt.size = 0.2, cols = col5)+
          coord_fixed()
ggsave("tissue.2-tsne_paired.pdf", device= "pdf", width = 20, 
               height = 16, units = "cm") 
        
save.image(file = object)
        
        
cell.type.order.2 <- c("epidermis_1",
"epidermis_2",
"epidermis_3",
"epidermis_4",
"epidermis_5",
"epidermis_6",
"nervous system_1",
"nervous system_2",
"nervous system_3",
"nervous system_4",
"nervous system_5",
"nervous system_6",
"nervous system_7",
"nervous system_8",
"nervous system_9",
"nervous system_10",
"nervous system_11",
"nervous system_12",
"nervous system_13",
"nervous system_14",
"nervous system_15",
"nervous system_16",
"nervous system_17",
"nervous system_18",
"nervous system_19",
"nervous system_20",
"nervous system_21",
"nervous system_22",
"endoderm_1",
"endoderm_2",
"endoderm_3",
"endoderm_4",
"mesenchyme_1",
"mesenchyme_2",
"mesenchyme_3",
"mesenchyme_4",
"mesenchyme_5",
"mesenchyme_6",
"mesenchyme_7",
"mesenchyme_8",
"trunk ventral cells",
"muscle",
"notochord",
"germ cells")
        
gfpvscrbn.combined$cell.type.2 <- factor(x= gfpvscrbn.combined$cell.type.2, levels = cell.type.order.2)

col6epi <- colorRampPalette(c("#A6CEE3" ,"#1F78B4"))(6)
col6NS <- colorRampPalette(c("#B2DF8A", "#33A02C"))(22)
col6endo <- colorRampPalette(c("#FB9A99", "#E31A1C"))(4)
col6mes <- colorRampPalette(c("#FDBF6F", "#FF7F00"))(8)
col6TVCmusc <- c("#CAB2D6", col4a[34])
col6germ<- col5germ
col6noto <- col5noto
              
col6 <- c(col6epi, col6NS, col6endo, col6mes, col6TVCmusc, col6noto, col6germ)

DimPlot(gfpvscrbn.combined, reduction = "tsne", 
                    group.by = "cell.type.2", label.size = 2, pt.size = 0.2,
                    label = FALSE, cols = col6)+
              coord_fixed()
ggsave("cell.types.2-tSNE_paired.pdf", device= "pdf", width = 32, 
                   height = 16, units = "cm") 
            
save.image(file = object)


# Cluster proportions between GFP and CRBN conditions

md2 <-subset(metadata, select = c("orig.ident", "seurat_clusters"))

is.data.table(md2)
md2 <- as.data.table(md2)

md2 <- md2 %>%group_by(seurat_clusters, orig.ident)

cluster.prop <- table(md2)
cluster.prop <- as.data.table(cluster.prop)
class(cluster.prop$seurat_clusters)
cluster.prop$seurat_clusters <- as.integer(cluster.prop$seurat_clusters)

cluster.prop1 <- count(md2, "orig.ident")
cluster.prop2 <- count(md2, "seurat_clusters")

#Proportion by cluster
for (x in 1:nrow(cluster.prop)) 
{
  a <- as.numeric(cluster.prop[x, "seurat_clusters"]) 
  print(a)
  b <- which(cluster.prop2$seurat_clusters %in% a)
  print(b)
  c <- as.numeric(cluster.prop2[b, n])
  print(c)
  cluster.prop[x, "prop"] <- as.numeric(cluster.prop[x, "N"]/c)
}
colnames(cluster.prop)[4]<-"prop.cluster"

#Proportion by origin
for (x in 1:nrow(cluster.prop)) 
{
  d <- as.character(cluster.prop[x, "orig.ident"]) 
  print(d)
  e <- which(cluster.prop1$orig.ident %in% d)
  print(e)
  f <- as.numeric(cluster.prop1[e, n])
  print(f)
  cluster.prop[x, "prop.orig"] <- as.numeric(cluster.prop[x, "N"]/f)
}

cluster.prop <- tidyft::left_join(cluster.prop,tissue.type2,
                                  by = c("seurat_clusters" = "cluster"))


cluster.prop$cell.type.2 <- factor(cluster.prop$cell.type.2, levels = cell.type.order.2)


Prop_graph1 <- ggplot(cluster.prop, 
                      aes(x = seurat_clusters, y = N, fill = orig.ident)) +
  geom_col()+
  ggtitle("Number of cells in each cluster")+
  scale_y_continuous(name ="Number of cells", expand = c(0,0))+
  scale_x_continuous(name = "Cluster #", labels = as.character(cluster.prop$seurat_clusters),
                     breaks = cluster.prop$seurat_clusters)+
  scale_fill_manual(values = rev(alpha(col1, 0.5)))


pdf(file="Cell.Number.by.Cluster.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
Prop_graph1+
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 8, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
dev.off()

#Proportion by cluster

Prop_graph2 <- ggplot(cluster.prop, 
                      aes(x = seurat_clusters, y = prop.cluster, fill = orig.ident)) +
  geom_col()+
  ggtitle("Proportion of cells in each cluster depending of the condition")+
  scale_y_continuous(name ="Proportion", expand = c(0,0))+
  scale_x_continuous(name = "Cluster #", labels = as.character(cluster.prop$seurat_clusters),
                     breaks = cluster.prop$seurat_clusters)+
  scale_fill_manual(values = rev(alpha(col1, 0.5)))


Prop_graph2c <- ggplot(cluster.prop, 
                       aes(x = cell.type.2, y = prop.cluster, fill = orig.ident)) + geom_col()+
  ggtitle("Cell Proportion in each cell type in function of the condition")+
  scale_y_continuous(name ="Proportion", expand = c(0,0))+
  scale_fill_manual(values = rev(alpha(col1, 0.5)))


pdf(file="Prop.norm to cluster size.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
Prop_graph2+
  #scale_x_continuous(labels = cluster.types)+
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
dev.off()



Prop_graph2c+
  #scale_x_continuous(labels = cluster.types)+
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 8, angle = 90, hjust = 1, vjust = 0.2),
        plot.title = element_text(hjust = 0.5))
ggsave("Ratio norm to cluster size2.pdf", device= "pdf", width = 20, 
       height = 10, units = "cm")

Prop_graph3 <- ggplot(cluster.prop, 
                      aes(x = orig.ident, y = prop.orig, fill = orig.ident)) +
  geom_bar(stat = "identity")+
  facet_grid(~seurat_clusters, switch = "x")+
  ggtitle("Cluster proportion in Tbx6>CRBN vs Tbx6>GFP")+
  scale_y_continuous(name ="Proportion", expand = c(0,0))+
  scale_x_discrete(name = "Cluster # by condition", labels = NULL,
                   breaks = NULL) +
  scale_fill_manual(name= "Condition",
                    values = rev(alpha(col1, 0.5)))

Prop_graph3b <- ggplot(cluster.prop, 
                      aes(x = orig.ident, y = prop.orig, fill = orig.ident)) +
  geom_bar(stat = "identity")+
  facet_grid(~cell.type.2, switch = "x")+
  ggtitle("Cluster proportion in Tbx6>CRBN vs Tbx6>GFP")+
  scale_y_continuous(name ="Proportion", expand = c(0,0))+
  scale_x_discrete(name = "Cell types", labels = NULL,
                   breaks = NULL) +
  scale_fill_manual(name= "Condition",
                    values = rev(alpha(col1, 0.5)))+
  theme_classic() +
  theme(strip.text.x =  element_text(color = "black", 
                                   size = 8, angle = -90, hjust = 0),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5))

pdf(file="Proportion Norm to origin.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
Prop_graph3+
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 6, angle = 45, hjust = 1),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()


pdf(file="Proportion Cell types Norm to origin.pdf", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
Prop_graph3b
dev.off()

save.image(object)


# Differentialy expressed genes in cluster expressing CFPSV40 
# (where CRBN is overexpressed)


DefaultAssay(gfpvscrbn.combined) <- "RNA"

pdf(file="CFPSV40.pdf", 
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches
VlnPlot(gfpvscrbn.combined, features = c("CFPSV40"), pt.size = 0,
        split.by = "orig.ident", group.by = "seurat_clusters", split.plot = TRUE, cols = alpha(col1, 0.5))
dev.off()




head(Idents(gfpvscrbn.combined))
Idents(gfpvscrbn.combined)<- "cell.type.2"

VlnPlot(gfpvscrbn.combined, features = c("CFPSV40"), pt.size = 0,
        split.by = "orig.ident", split.plot = TRUE, cols = alpha(col1, 0.5))+
  ggtitle("GFP reporter")+
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 8, angle = 90, hjust = 1, vjust = 0.2),
        plot.title = element_text(hjust = 0.5))
ggsave("GFP expression level2.pdf", device= "pdf", width = 20, 
       height = 8, units = "cm")

p1 <- VlnPlot(gfpvscrbn.combined, features = c("CFPSV40"), pt.size = 0,
        split.by = "orig.ident", split.plot = TRUE, cols = alpha(col1, 0.5))+
  ggtitle("GFP reporter")+
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 8, angle = 90, hjust = 1, vjust = 0.2),
        plot.title = element_text(hjust = 0.5))+
  coord_fixed()

p1$layers[[1]]$aes_params$size = 0.2 #Change thickness of plot



pdf(file="GFP expression level3.pdf", 
    width = 10, # The width of the plot in inches
    height = 6) # The height of the plot in inches
p1
dev.off()


save.image(object)


# tested clusters: 5 (mesenchyme), 10 (mesenchyme), 11 (endoderm), 13 (muscle), 18 (mesenchyme)
# 25 (mesenchyme), 39 (TVCs), 41 (mesenchyme)

# For performing differential expression, we switch back to the original data (RNA assay)

DefaultAssay(gfpvscrbn.combined) <- "RNA"  

# Normalized and scale the RNA assay of the combined data set 
#(the RNA assay is use to find DE genes, the scale data are used to plot heatmap in seurat,
# not necessary if used complex heatmap)

all.genes <- rownames(gfpvscrbn.combined)
gfpvscrbn.combined <- ScaleData(gfpvscrbn.combined, features = all.genes)

head(gfpvscrbn.combined@meta.data)
head(Idents(gfpvscrbn.combined))

#Create metadata with combined cluster ID and electroporation conditions

gfpvscrbn.combined@meta.data$GroupCluster <- paste(gfpvscrbn.combined@meta.data$seurat_clusters,
                                                   gfpvscrbn.combined@meta.data$orig.ident, sep = "_")


Idents(gfpvscrbn.combined) <- "GroupCluster" #modify identity

g <- c(5, 10, 11, 13, 18, 25 , 39, 41) #select the cluster to look for DE genes

cluster.de <- lapply(g, function(cl){
  
  cl.de<- FindMarkers(gfpvscrbn.combined, ident.1 = paste(cl,"Tbx6_CRBN", sep ="_"),
                      ident.2 = paste(cl,"Tbx6_GFP", sep ="_"), 
                      min.pct = 0.25,
                      logfc.threshold = log2(1.5))
  cl.de$cluster <- cl
  cl.de$gene <- rownames(cl.de)
  return(cl.de)
})

cluster.de <- bind_rows(cluster.de)
is.data.table(cluster.de)

cluster.de <- tidyft::left_join(cluster.de,human.homo,
                                by = c("gene" = "KH.model.ID"))
write.table(cluster.de, file="CRBN.OE.de_filteredMatrices_CCA44.txt",
            quote=F, sep="\t", col.names=NA)

save.image(file = object)

#If use 0.5 log 2 fold change

cluster.de.0.5 <- lapply(g, function(cl){
  
  cl.de<- FindMarkers(gfpvscrbn.combined, ident.1 = paste(cl,"Tbx6_CRBN", sep ="_"),
                      ident.2 = paste(cl,"Tbx6_GFP", sep ="_"), 
                      min.pct = 0.25,
                      logfc.threshold = 0.5)
  cl.de$cluster <- cl
  cl.de$gene <- rownames(cl.de)
  return(cl.de)
})

cluster.de.0.5 <- bind_rows(cluster.de.0.5)
is.data.table(cluster.de.0.5)

cluster.de.0.5 <- tidyft::left_join(cluster.de.0.5,human.homo,
                                by = c("gene" = "KH.model.ID"))
write.table(cluster.de.0.5, file="CRBN.OE.de_filteredMatrices_CCA44_0.5.txt",
            quote=F, sep="\t", col.names=NA)

save.image(file = object)

# Add cell type to degs:  5 (mesenchyme), 10 (mesenchyme), 11 (endoderm), 13 (muscle), 18 (mesenchyme)
# 25 (mesenchyme), 39 (TVCs), 41 (mesenchyme)

cluster.de.0.5_2 <- cluster.de.0.5

for (x in 1:nrow(cluster.de.0.5))
{
  cluster.de.0.5[x, "cell.type"] <- if (cluster.de.0.5[x, "cluster"] == "5")
  {"mesenchyme_1"} 
  else 
  { if (cluster.de.0.5[x, "cluster"] == "10")
  {"mesenchyme_2"}
    else
    { if (cluster.de.0.5[x, "cluster"] == "11")
    {"endoderm_2"}
      else { if (cluster.de.0.5[x, "cluster"] == "13")
      {"muscle"}
        else { if (cluster.de.0.5[x, "cluster"] == "18")
        {"mesenchyme_4"}
          else {if(cluster.de.0.5[x, "cluster"] == "25")
          {"mesenchyme_6"}
            else {if (cluster.de.0.5[x, "cluster"] == "39")
            {"trunk ventral cells"}
              else {if (cluster.de.0.5[x, "cluster"] == "41")
              {"mesenchyme_8"}
              }
            }
          }
          }
        }
       }
    }
  }


write.table(cluster.de.0.5, file="CRBN.OE.de_filteredMatrices_CCA44_0.5_annotated.txt",
            quote=F, sep="\t", col.names=NA)

# subset deg for only p val adjust < 0.5

cluster.de.sign.fc0.5 <- cluster.de.0.5[ which(cluster.de.0.5$p_val_adj < 0.05), ]

cluster.de.sign <- cluster.de[ which(cluster.de$p_val_adj < 0.05), ]
  
# Heatmap_muscle-cluster13

top.de.cluster13 <- cluster.de.sign.fc0.5[ which(cluster.de.sign.fc0.5$cluster == 13), ]

is.data.table(top.de.cluster13)
top.de.cluster13 <- as.data.table(top.de.cluster13)

top.de.cluster13<-top.de.cluster13%>%mutate(Human.homolog=coalesce(Human.homolog,gene))

name.top.de.cluster13 <- structure(as.character(top.de.cluster13$Human.homolog),
                                   names = as.character(top.de.cluster13$gene))

head(name.top.de.cluster13)

cluster.cells.13 <- subset(metadata, seurat_clusters == 13,
                           select=c(orig.ident, seurat_clusters))

cluster.cells.13 <- rownames_to_column(cluster.cells.13)

head(cluster.cells.13)

mat13 <- gfpvscrbn.combined[["RNA"]]@data[top.de.cluster13$gene,
                                          cluster.cells.13$rowname]
mat13 <- as.matrix(mat13)
mat13<- t(scale(t(mat13)))
quantile(mat13, c(0.1, 0.95))

anno<- data.frame(cluster.cells.13$orig.ident,cluster.cells.13$seurat_clusters)

colnames(anno)<-c("condition", "cluster")

column_ha = HeatmapAnnotation(df = anno,
                              col = list("condition" = c("Tbx6_GFP" = alpha(col1[1],0.5),
                                                         "Tbx6_CRBN" = alpha(col1[2],0.5)),
                                         "cluster" = c("13" = "lightblue1")))


#Heatmap for muscle cell (cluster13) with genes order based on fold changed and split samples

top.de.cluster13.name <- subset(top.de.cluster13, select = c("avg_log2FC","gene", "Human.homolog"))

top.de.cluster13.name <- top.de.cluster13.name[order(-top.de.cluster13.name$avg_log2FC)]
order3 <- top.de.cluster13.name$gene

cluster.cells.13 <- cluster.cells.13[c(1:3)]
cluster.cells.13.2 <- column_to_rownames(cluster.cells.13, var = "rowname")
head(cluster.cells.13.2)

display.brewer.pal(11, "RdBu")
brewer.pal(11, "RdBu")

col_fun = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

pdf(file="Heatmap_anno-split-orderfoldchanged_cluster13_fc1.4.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
Heatmap(mat13, name = "Rel Expr", 
        row_labels = name.top.de.cluster13[rownames(mat13)],
        row_names_gp = gpar(fontsize = 7),
        column_split = factor(cluster.cells.13.2$orig.ident, levels = c("Tbx6_GFP", "Tbx6_CRBN")),
        column_title = c("Tbx6>GFP", "Tbx6>Crbn"),
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        row_order = order3,
        show_column_names = FALSE,
        bottom_annotation = column_ha,
        col = col_fun,
        width = ncol(mat13)*unit(0.2, "mm"),
        height = nrow(mat13)*unit(3, "mm"))
dev.off() 

pdf(file="Heatmap_split-orderfoldchanged_cluster13_fc1.4.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
Heatmap(mat13, name = "Rel Expr", 
        #row_labels = name.top.de.cluster13[rownames(mat13)],
        row_names_gp = gpar(fontsize = 7),
        column_split = factor(cluster.cells.13.2$orig.ident, levels = c("Tbx6_GFP", "Tbx6_CRBN")),
        column_title = c("Tbx6>GFP", "Tbx6>Crbn"),
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        row_order = order3,
        show_column_names = FALSE,
        bottom_annotation = column_ha,
        col = col_fun,
        width = ncol(mat13)*unit(0.2, "mm"),
        height = nrow(mat13)*unit(3, "mm"))
dev.off() 

save.image(file=object)

#Heatmap cluster 13 with 1.5 fold change 

top.1.5.de.cluster13 <- cluster.de.sign[ which(cluster.de.sign$cluster == 13), ]
is.data.table(top.1.5.de.cluster13)
top.1.5.de.cluster13 <- as.data.table(top.1.5.de.cluster13)
top.1.5.de.cluster13<-top.1.5.de.cluster13%>%mutate(Human.homolog=coalesce(Human.homolog,gene))
name.top.1.5.de.cluster13 <- structure(as.character(top.1.5.de.cluster13$Human.homolog),
                                   names = as.character(top.1.5.de.cluster13$gene))

head(name.top.1.5.de.cluster13)

mat13.fc1.5 <- gfpvscrbn.combined[["RNA"]]@data[top.1.5.de.cluster13$gene,
                                          cluster.cells.13$rowname]
mat13.fc1.5 <- as.matrix(mat13.fc1.5)
mat13.fc1.5<- t(scale(t(mat13.fc1.5)))
quantile(mat13.fc1.5, c(0.1, 0.95))

#Heatmap for muscle cell (cluster13) with genes order based on fold changed and split samples_fold change 1.5

top.1.5.de.cluster13.name <- subset(top.1.5.de.cluster13, select = c("avg_log2FC","gene", "Human.homolog"))

top.1.5.de.cluster13.name <- top.1.5.de.cluster13.name[order(-top.1.5.de.cluster13.name$avg_log2FC)]
order3.fc1.5 <- top.1.5.de.cluster13.name$gene

pdf(file="Heatmap_anno-split-orderfoldchanged_cluster13_fc1.5.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
Heatmap(mat13.fc1.5, name = "Rel Expr", 
        row_labels = name.top.1.5.de.cluster13[rownames(mat13.fc1.5)],
        row_names_gp = gpar(fontsize = 7),
        column_split = factor(cluster.cells.13.2$orig.ident, levels = c("Tbx6_GFP", "Tbx6_CRBN")),
        column_title = c("Tbx6>GFP", "Tbx6>Crbn"),
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        row_order = order3.fc1.5,
        show_column_names = FALSE,
        bottom_annotation = column_ha,
        col = col_fun,
        width = ncol(mat13)*unit(0.2, "mm"),
        height = nrow(mat13)*unit(3, "mm"))
dev.off() 

pdf(file="Heatmap_split-orderfoldchanged_cluster13_fc1.5.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
Heatmap(mat13.fc1.5, name = "Rel Expr", 
        #row_labels = name.top.1.5.de.cluster13[rownames(mat13)],
        row_names_gp = gpar(fontsize = 7),
        column_split = factor(cluster.cells.13.2$orig.ident, levels = c("Tbx6_GFP", "Tbx6_CRBN")),
        column_title = c("Tbx6>GFP", "Tbx6>Crbn"),
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        row_order = order3.fc1.5,
        show_column_names = FALSE,
        bottom_annotation = column_ha,
        col = col_fun,
        width = ncol(mat13)*unit(0.2, "mm"),
        height = nrow(mat13)*unit(3, "mm"))
dev.off() 

save.image(file=object)

## Heatmap for cluster 10 (mesenchyme)_1.5 fold change

top.de.cluster10 <- cluster.de.sign[ which(cluster.de.sign$cluster == 10), ]

is.data.table(top.de.cluster10)
top.de.cluster10 <- as.data.table(top.de.cluster10)

top.de.cluster10<-top.de.cluster10%>%mutate(Human.homolog=coalesce(Human.homolog,gene))

name.top.de.cluster10 <- structure(as.character(top.de.cluster10$Human.homolog),
                                   names = as.character(top.de.cluster10$gene))

head(name.top.de.cluster10)

cluster.cells.10 <- subset(metadata, seurat_clusters == 10,
                           select=c(orig.ident, seurat_clusters))

cluster.cells.10 <- rownames_to_column(cluster.cells.10)


mat10 <- gfpvscrbn.combined[["RNA"]]@data[top.de.cluster10$gene,
                                          cluster.cells.10$rowname]
mat10 <- as.matrix(mat10)
mat10<- t(scale(t(mat10)))
quantile(mat10, c(0.1, 0.95))

anno10<- data.frame(cluster.cells.10$orig.ident,cluster.cells.10$seurat_clusters)

colnames(anno10)<-c("condition", "cluster")

column_ha10 = HeatmapAnnotation(df = anno10,
                              col = list("condition" = c("Tbx6_GFP" = alpha(col1[1],0.5),
                                                         "Tbx6_CRBN" = alpha(col1[2],0.5)),
                                         "cluster" = c("10" = "violet")))


#Heatmap for cluster10 with genes order based on fold changed and split samples

top.de.cluster10.name <- subset(top.de.cluster10, select = c("avg_log2FC","gene", "Human.homolog"))

top.de.cluster10.name <- top.de.cluster10.name[order(-top.de.cluster10.name$avg_log2FC)]

order4 <- top.de.cluster10.name$gene

cluster.cells.10.2 <- column_to_rownames(cluster.cells.10, var = "rowname")

head(cluster.cells.10.2)

#col_fun2 = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

pdf(file="Heatmap_anno-split-orderfoldchanged_cluster10_fc1.5.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
Heatmap(mat10, name = "Rel Expr", 
        row_labels = name.top.de.cluster10[rownames(mat10)],
        row_names_gp = gpar(fontsize = 7),
        column_split = factor(cluster.cells.10.2$orig.ident, levels = c("Tbx6_GFP", "Tbx6_CRBN")),
        column_title = c("Tbx6>GFP", "Tbx6>Crbn"),
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        row_order = order4,
        show_column_names = FALSE,
        bottom_annotation = column_ha10,
        col = col_fun,
        width = ncol(mat10)*unit(0.2, "mm"),
        height = nrow(mat10)*unit(3, "mm"))
dev.off() 

pdf(file="Heatmap_split-orderfoldchanged_cluster10_fc1.5.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
Heatmap(mat10, name = "Rel Expr", 
        #row_labels = name.top.de.cluster10[rownames(mat10)],
        row_names_gp = gpar(fontsize = 7),
        column_split = factor(cluster.cells.10.2$orig.ident, levels = c("Tbx6_GFP", "Tbx6_CRBN")),
        column_title = c("Tbx6>GFP", "Tbx6>Crbn"),
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        row_order = order4,
        show_column_names = FALSE,
        bottom_annotation = column_ha10,
        col = col_fun,
        width = ncol(mat10)*unit(0.2, "mm"),
        height = nrow(mat10)*unit(3, "mm"))
dev.off() 

save.image(file=object)


## Cluster 18 (mesenchyme)_1.5 fold change

top.de.cluster18 <- cluster.de.sign[ which(cluster.de.sign$cluster == 18), ]

is.data.table(top.de.cluster18)
top.de.cluster18 <- as.data.table(top.de.cluster18)

top.de.cluster18<-top.de.cluster18%>%mutate(Human.homolog=coalesce(Human.homolog,gene))

name.top.de.cluster18 <- structure(as.character(top.de.cluster18$Human.homolog),
                                   names = as.character(top.de.cluster18$gene))

head(name.top.de.cluster18)

cluster.cells.18 <- subset(metadata, seurat_clusters == 18,
                           select=c(orig.ident, seurat_clusters))

cluster.cells.18 <- rownames_to_column(cluster.cells.18)

mat18 <- gfpvscrbn.combined[["RNA"]]@data[top.de.cluster18$gene,
                                          cluster.cells.18$rowname]
mat18 <- as.matrix(mat18)
mat18<- t(scale(t(mat18)))
quantile(mat18, c(0.1, 0.95))

anno18<- data.frame(cluster.cells.18$orig.ident,cluster.cells.18$seurat_clusters)

colnames(anno18)<-c("condition", "cluster")

column_ha18 = HeatmapAnnotation(df = anno18,
                                col = list("condition" = c("Tbx6_GFP" = alpha(col1[1],0.5),
                                                           "Tbx6_CRBN" = alpha(col1[2],0.5)),
                                           "cluster" = c("18" = "violet")))


#Heatmap for cluster18 with genes order based on fold changed and split samples

top.de.cluster18.name <- subset(top.de.cluster18, select = c("avg_log2FC","gene", "Human.homolog"))

top.de.cluster18.name <- top.de.cluster18.name[order(-top.de.cluster18.name$avg_log2FC)]

order5 <- top.de.cluster18.name$gene

cluster.cells.18.2 <- column_to_rownames(cluster.cells.18, var = "rowname")

head(cluster.cells.18.2)

#col_fun2 = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

pdf(file="Heatmap_anno-split-orderfoldchanged_cluster18_fc1.5.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
Heatmap(mat18, name = "Rel Expr", 
        row_labels = name.top.de.cluster18[rownames(mat18)],
        row_names_gp = gpar(fontsize = 7),
        column_split = factor(cluster.cells.18.2$orig.ident, levels = c("Tbx6_GFP", "Tbx6_CRBN")),
        column_title = c("Tbx6>GFP", "Tbx6>Crbn"),
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        row_order = order5,
        show_column_names = FALSE,
        bottom_annotation = column_ha18,
        col = col_fun,
        width = ncol(mat18)*unit(0.2, "mm"),
        height = nrow(mat18)*unit(3, "mm"))
dev.off() 

pdf(file="Heatmap_split-orderfoldchanged_cluster18_fc1.5.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
Heatmap(mat18, name = "Rel Expr", 
        #row_labels = name.top.de.cluster18[rownames(mat18)],
        row_names_gp = gpar(fontsize = 7),
        column_split = factor(cluster.cells.18.2$orig.ident, levels = c("Tbx6_GFP", "Tbx6_CRBN")),
        column_title = c("Tbx6>GFP", "Tbx6>Crbn"),
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        row_order = order5,
        show_column_names = FALSE,
        bottom_annotation = column_ha18,
        col = col_fun,
        width = ncol(mat18)*unit(0.2, "mm"),
        height = nrow(mat18)*unit(3, "mm"))
dev.off() 

save.image(file=object)

## Heatmap for cluster 25 _fold change 1.5

top.de.cluster25 <- cluster.de.sign[ which(cluster.de.sign$cluster == 25), ]

is.data.table(top.de.cluster25)
top.de.cluster25 <- as.data.table(top.de.cluster25)

top.de.cluster25<-top.de.cluster25%>%mutate(Human.homolog=coalesce(Human.homolog,gene))

name.top.de.cluster25 <- structure(as.character(top.de.cluster25$Human.homolog),
                                   names = as.character(top.de.cluster25$gene))

head(name.top.de.cluster25)

cluster.cells.25 <- subset(metadata, seurat_clusters == 25,
                           select=c(orig.ident, seurat_clusters))

cluster.cells.25 <- rownames_to_column(cluster.cells.25)

mat25 <- gfpvscrbn.combined[["RNA"]]@data[top.de.cluster25$gene,
                                          cluster.cells.25$rowname]
mat25 <- as.matrix(mat25)
mat25<- t(scale(t(mat25)))
quantile(mat25, c(0.1, 0.95))

anno25<- data.frame(cluster.cells.25$orig.ident,cluster.cells.25$seurat_clusters)

colnames(anno25)<-c("condition", "cluster")

column_ha25 = HeatmapAnnotation(df = anno25,
                                col = list("condition" = c("Tbx6_GFP" = alpha(col1[1],0.5),
                                                           "Tbx6_CRBN" = alpha(col1[2],0.5)),
                                           "cluster" = c("25" = "violet")))


#Heatmap for cluster25 with genes order based on fold changed and split samples

top.de.cluster25.name <- subset(top.de.cluster25, select = c("avg_log2FC","gene", "Human.homolog"))

top.de.cluster25.name <- top.de.cluster25.name[order(-top.de.cluster25.name$avg_log2FC)]

order6 <- top.de.cluster25.name$gene

cluster.cells.25.2 <- column_to_rownames(cluster.cells.25, var = "rowname")

head(cluster.cells.25.2)

#col_fun2 = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

pdf(file="Heatmap_anno-split-orderfoldchanged_cluster25_fc1.5.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
Heatmap(mat25, name = "Rel Expr", 
        row_labels = name.top.de.cluster25[rownames(mat25)],
        row_names_gp = gpar(fontsize = 7),
        column_split = factor(cluster.cells.25.2$orig.ident, levels = c("Tbx6_GFP", "Tbx6_CRBN")),
        column_title = c("Tbx6>GFP", "Tbx6>Crbn"),
        cluster_column_slices = FALSE,
        show_column_dend = FALSE,
        row_order = order6,
        show_column_names = FALSE,
        bottom_annotation = column_ha25,
        col = col_fun,
        width = ncol(mat25)*unit(0.2, "mm"),
        height = nrow(mat25)*unit(3, "mm"))
dev.off() 

save.image(file=object)

#Check if notochord also show signs of more immature cells:

Idents(gfpvscrbn.combined) <- "GroupCluster" 
head(Idents(gfpvscrbn.combined))

noto_deg<- FindMarkers(gfpvscrbn.combined, ident.1 = "9_Tbx6_CRBN",
                       ident.2 = "9_Tbx6_GFP", 
                       min.pct = 0.25,
                       logfc.threshold = 0.5)

noto_deg <- rownames_to_column(noto_deg)
colnames(noto_deg)[1] <- "gene"

noto_deg <-tidyft::left_join(noto_deg,human.homo,
                             by = c("gene" = "KH.model.ID"))
save.image(file = object)


#Check number of cells

samplename = gfpvscrbn.combined@meta.data$orig.ident
table(samplename)

save.image(file = object)

##Check expression level of Myod in muscle cells

Idents(gfpvscrbn.combined) <- "cell.type.2"

VlnPlot(gfpvscrbn.combined, features = c("KH2012:KH.C14.307"), pt.size = 0.5,
        split.by = "orig.ident", split.plot = TRUE, idents = "muscle", cols = alpha(col1, 0.5))+
  # ggtitle("GFP reporter")+
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 8, angle = 90, hjust = 1, vjust = 0.2),
        plot.title = element_text(hjust = 0.5))
ggsave("Violin Plot Myod.pdf", device= "pdf", width = 12, 
       height = 8, units = "cm")


matMyod <- gfpvscrbn.combined[["RNA"]]@data[c("KH2012:KH.C14.307"),
                                                 cluster.cells.13$rowname]

matMyod <- as.matrix(matMyod)
matMyod <- as.data.frame(matMyod)

matMyod <- rownames_to_column(matMyod)
matMyod[c('cell', 'orig.ident')] <- str_split_fixed(matMyod$rowname, '-', 2)


matMyod = subset(matMyod, select = -c(cell) )
colnames(matMyod)[1] <- "Cell"
colnames(matMyod)[2] <- "ExpressionLevel"


matMyod$Gene <- "Myod"

matMyod$orig.ident <- factor(matMyod$orig.ident, levels = c("Tbx6_GFP", "Tbx6_CRBN"))

ggplot(matMyod, aes(x = orig.ident, y = ExpressionLevel, color = orig.ident)) + 
  geom_boxplot(aes(fill = orig.ident), outlier.shape = NA,
               color = alpha(col1, 0.8))+
  scale_fill_manual(values =  rep("white", each = 2))+
  geom_point(position=position_jitterdodge(), size = 0.2)+
  scale_color_manual(values = rep("black", each = 2))+
  
  ggtitle("Myod expression level")+
  scale_y_continuous(limits = c(0, 2.2) )+
  labs( y = "Log Normalized Expression Level" , x= "") + 
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 12),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/2)

ggsave("Box+dot Plot Myod.pdf", device= "pdf",
       width = 12, height = 8, units = "cm")


save.image(object)

##Statistical test for Myod expression 

## test if normal distribution
shapiro.test(matMyod$ExpressionLevel)

#When the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. 
#In other words, we can assume the normality.
#It is say to be normal but it could be not enough datapoint to test

t.test(ExpressionLevel ~ orig.ident, data = matMyod, alternative = "two.sided")

## p-value <0.05, the distribution is significantly different from normal
## Mann-Withney U test (Wilcox test)

wilcox.test(ExpressionLevel ~ orig.ident, data = matMyod, exact = FALSE)


save.image(file=object)

sessionInfo()
