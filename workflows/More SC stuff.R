library(Seurat)
library(SeuratData)
library(patchwork)

LoadData("bmcite")

bmcite

cell.dat <- spectre()
cell.dat

cell.dat@meta <- data.table('CellID' = colnames(bmcite))
cell.dat@data$RNA <- Matrix::t(bmcite@assays$RNA@counts)
cell.dat@data$ADT <- as.data.table(Matrix::t(bmcite@assays$ADT@counts))

cell.dat@data$RNA
cell.dat


###

is.mito <- grepl("^MT-", colnames(cell.dat@data$RNA))
is.mito <- colnames(cell.dat@data$RNA)[is.mito]
is.mito

qcstats <- perCellQCMetrics(Matrix::t(cell.dat@data$RNA), subsets=list(Mito=is.mito))
qcstats

filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
filtered

!filtered$discard # rows for exclusion

### Filtered

cell.dat <- sc.row.filter(cell.dat, !filtered$discard)
cell.dat

### Normalisation

cell.dat <- sc.do.lognorm(cell.dat, 'RNA')
cell.dat

### Feat selection

dec <- scran::modelGeneVar(Matrix::t(cell.dat@data$RNA_logNorm))
hvg <- scran::getTopHVGs(dec, prop=0.1)
hvg

cell.dat@info$hvg <- hvg
cell.dat@info

### PCA

cell.dat <- sc.run.pca(cell.dat, 'RNA_logNorm', use.cols = cell.dat@info$hvg)
cell.dat

### Run UMAP

cell.dat <- sc.run.umap(cell.dat, 'PCA', slot = 'analysis')
cell.dat

cell.dat@






# 
# 
# make.colour.plot(dat, 
#                  
#                  c('UMAP_1', 'UMAP'),
#                  c('UMAP_2', 'UMAP'),
#                  c('CD3', 'asinh'),
#                  
#                  x.slot = NULL,
#                  y.slot = NULL,
#                  c.slot = NULL
#                  
#                  )
# 
# # dat
# 
# ###
#     # CD3
#     # CD1c
#     # CD4
# 
# ###
#     # asinh
#     # asinh
#     # asinh
# 
# ###
#     # data
#     # data
#     # data



