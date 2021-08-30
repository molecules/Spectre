#' sc.row.filter
#' 
#' @export

sc.row.filter <- function(dat, rows){
  
  nrow.start <- nrow(dat@data[[1]])
  
  res <- spectre()
  res@meta <- dat@meta[rows,]
  
  for(i in names(dat@data)){
    res@data[[i]] <- dat@data[[i]][rows,]
  }
  
  for(i in names(dat@analysis)){
    res@analysis[[i]] <- dat@analysis[[i]][rows,]
  }
  
  nrow.end <- nrow(res@data[[1]])
  
  message(' -- Rows before: ', nrow.start)
  message(' -- Rows after: ', nrow.end)
  
  return(res)
}



#' sc.do.lognorm
#' 
#' @export

sc.do.lognorm <- function(dat, target, use.cols){
  
  require('Seurat')
  require('data.table')
  
  res <- Seurat::LogNormalize(Matrix::t(dat@data[[target]]))
  new.name <- paste0(target, '_logNorm')
  dat@data[[new.name]] <- Matrix::t(res)
  
  return(dat)
}




#' sc.run.pca
#' 
#' @export

sc.run.pca <- function(dat, target, use.cols = NULL, npcs = 30, set.seed = 42){
  
  require('scater')
  require('data.table')
  
  # dat <- new
  # name <- 'RNA_logNorm'
  # use.cols <- hvg
  # npcs = 30
  # set.seed = 42
  
  set.seed(set.seed)
  
  if(is.null(use.cols)){
    use.cols <- colnames(dat@data[[target]])
  }
  
  res <- scater::calculatePCA(Matrix::t(dat@data[[target]]), ncomponents=npcs, subset_row=use.cols)
  res <- as.data.table(res)
  
  dat@analysis[['PCA']] <- res
  return(dat)
}


#' sc.run.umap
#' 
#' @export

sc.run.umap <- function(dat, target, use.cols = NULL, slot = 'data'){
  
  # dat <- cell.dat
  # name <- 'PCA'
  # use.cols <- NULL
  # slot <- 'analysis'
  
  which(target == names(dat@data))
  which(target == names(dat@analysis))
  
  if(slot == 'data'){
    dat@data[[target]]
    if(is.null(use.cols)){
      use.cols <- names(dat@data[[target]])
    }
    res <- Spectre::run.umap(dat@data[[target]], use.cols)
  }
  
  if(slot == 'analysis'){
    dat@analysis[[target]]
    if(is.null(use.cols)){
      use.cols <- names(dat@analysis[[target]])
    }
    
    res <- Spectre::run.umap(dat@analysis[[target]], use.cols)
    x <- length(names(res))
    res <- res[,c(x-1, x), with = FALSE]
  }
  
  dat@analysis[['UMAP']] <- res
  return(dat)
}



