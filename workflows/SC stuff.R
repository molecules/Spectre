library(Seurat)
library(SeuratData)
library(patchwork)

# 
# # install dataset
# InstallData("ifnb")
# # load dataset
# LoadData("ifnb")
# 
# # split the dataset into a list of two seurat objects (stim and CTRL)
# ifnb.list <- SplitObject(ifnb, split.by = "stim")
# 
# 
# ifnb.list





# 
# # normalize and identify variable features for each dataset independently
# ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# 
# # select features that are repeatedly variable across datasets for integration
# features <- SelectIntegrationFeatures(object.list = ifnb.list)
# 



library('Spectre')


library(dplyr)
library(Seurat)
library(patchwork)
# 
# # Load the PBMC dataset
# pbmc.data <- Read10X(data.dir = "~/OneDrive - The University of Sydney (Staff)/Library/Data/filtered_gene_bc_matrices/hg19/")
# pbmc.data
# 
# # Initialize the Seurat object with the raw (non-normalized data).
# pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# pbmc




# 

cell.dat <- Spectre::demo.batches.1

cell.dat.a <- do.filter(cell.dat, 'Batch', 'A')
cell.dat.a <- do.subsample(cell.dat.a, 1000)

cell.dat.b <- do.filter(cell.dat, 'Batch', 'B')
cell.dat.b <- do.subsample(cell.dat.b, 1000)

# plot(cell.dat.b$`DL800 Ly6G_asinh`)
# cell.dat.b <- cell.dat.b[cell.dat.b[['DL800 Ly6G_asinh']] < 1.5,]






as.matrix(names(cell.dat.a))

meta.a <- cell.dat.a[, names(cell.dat.a)[c(19:21)], with = FALSE]
meta.a

as.matrix(names(cell.dat.a))

counts.a <- cell.dat.a[, names(cell.dat.a)[c(11:18)], with = FALSE]
counts.a <- as.matrix(counts.a)
counts.a <- Matrix::Matrix(counts.a, sparse = TRUE)
rownames(counts.a) <- paste0('CellA-', c(1:nrow(counts.a)))
counts.a <- t(counts.a)




as.matrix(names(cell.dat.b))

meta.b <- cell.dat.b[, names(cell.dat.b)[c(19:21)], with = FALSE]
meta.b

as.matrix(names(cell.dat.b))

counts.b <- cell.dat.b[, names(cell.dat.b)[c(11:18)], with = FALSE]
counts.b <- as.matrix(counts.b)
counts.b <- Matrix::Matrix(counts.b, sparse = TRUE)
rownames(counts.b) <- paste0('CellB-', c(1:nrow(counts.b)))
counts.b <- t(counts.b)


# counts <- dat[,names(dat)[c(11:19)], with = FALSE]
# counts <- as.matrix(counts)
# counts <- Matrix::Matrix(counts, sparse = TRUE)
# counts
# 
# rownames(counts) <- paste0('Cell-', c(1:nrow(counts)))
# counts <- t(counts)


srt.a <- CreateSeuratObject(counts = counts.a, assay = 'cyto')
srt.a <- FindVariableFeatures(srt.a, selection.method = "vst", nfeatures = 8)
srt.a

srt.b <- CreateSeuratObject(counts = counts.b, assay = 'cyto')
srt.b <- FindVariableFeatures(srt.b, selection.method = "vst", nfeatures = 8)
srt.b

features <- SelectIntegrationFeatures(object.list = list(srt.a, srt.b))
features


#####

# srt.a <- ScaleData(srt.a, features = features, verbose = FALSE, assay = 'cyto')
# srt.b <- ScaleData(srt.b, features = features, verbose = FALSE, assay = 'cyto')
# 
# srt.a@assays$cyto@scale.data
# 
# srt.a <- RunPCA(srt.a, features = features, verbose = FALSE, assay = 'cyto')
# srt.b <- RunPCA(srt.b, features = features, verbose = FALSE, assay = 'cyto')
# 
# srt.a@reductions$pca

Seurat::RunCCA()


immune.anchors <- FindIntegrationAnchors(object.list = list(srt.a, srt.b), 
                                         anchor.features = features, 
                                         dims = 1:7, 
                                         reduction = 'cca')#, reduction = "rpca")
immune.anchors

immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:7)
DefaultAssay(immune.combined) <- "integrated"


setwd('~/Desktop')


x <- immune.combined@assays$integrated@data

x <- t(x)
x <- as.data.table(x)
x <- cbind(x, rbind(meta.a, meta.b))
x

as.matrix(names(x))

x <- run.fitsne(x, names(x)[c(1:8)])
setwd('~/Desktop')
make.colour.plot(x, 'FItSNE_X', 'FItSNE_Y', 'Batch', filename = 'Integrated.png')
make.multi.plot(x, 'FItSNE_X', 'FItSNE_Y', names(x)[c(1:8)], figure.title = 'Integrated multi')

  
  
  
  
y <- immune.combined@assays$cyto@data

y <- t(y)
y <- as.data.table(y)
y <- cbind(y, rbind(meta.a, meta.b))
y

as.matrix(names(y))

y <- run.fitsne(y, names(y)[c(1:8)])
make.colour.plot(y, 'FItSNE_X', 'FItSNE_Y', 'Batch', filename = 'Raw.png')
make.multi.plot(y, 'FItSNE_X', 'FItSNE_Y', names(y)[c(1:8)], figure.title = 'Raw multi')









#run.cca

run.cca # with mnn
run.rcpa # with mnn
run.harmony # no mnn


run.mnn-adjust






do.split <- function(dat,
                     by,
                     type = 'meta'){
  
  ###
  
  res.list <- list()
  
  ###
  
  if(type == 'meta'){
    lst <- unique(dat@meta[[by]])
  }
  
  ###
  
  rws.a <- dat@meta[['Batch']] == 'A'
  rws.b <- dat@meta[['Batch']] == 'B'
  
  dat.a <- dat@data$asinh[rws.a,]
  dat.b <- dat@data$asinh[rws.b,]
  
  
  ###
  
  
}

dat.list <- do.split(cell.dat, 'Batch')









