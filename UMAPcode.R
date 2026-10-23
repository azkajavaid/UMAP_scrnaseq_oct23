library(Matrix)
library(Seurat)

control <- as.matrix(readMM("GSM2560248_2.1.mtx")); 
barcodes_ctrl <- read.table("GSM2560248_barcodes.tsv");
colnames(control) <- paste("ctrl_", barcodes_ctrl$V1, sep=""); 
rownames(control) <- read.table("GSE96583_batch2.genes.tsv")$V2;

stim <- as.matrix(readMM("GSM2560249_2.2.mtx")); 
barcodes_stim <- read.table("GSM2560249_barcodes.tsv");
colnames(stim) <- paste("stim_", barcodes_stim$V1, sep=""); 
rownames(stim) <- read.table("GSE96583_batch2.genes.tsv")$V2

tot <- t(rbind(t(control), t(stim)))
rownames(tot) <- make.unique(rownames(tot))
attr <- read.delim("GSE96583_batch2.total.tsne.df.tsv", header = TRUE, row.names = 1, check.names = FALSE)
identical(rownames(attr), c(barcodes_ctrl$V1, barcodes_stim$V1))
seurat_obj <- CreateSeuratObject(counts = tot)
seurat_obj <- NormalizeData(seurat_obj)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap")

umap_coords <- as.data.frame(seurat_obj[["umap"]]@cell.embeddings)
do.call(data.frame, strsplit(colnames(seurat_obj), "_")[1,])
cbind(umap_coords, seurat_obj$seurat_clusters)

