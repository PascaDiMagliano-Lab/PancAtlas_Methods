generate_cell_density <- function(seurat, col_name, cluster_name, weights, ...){
  coords <- seurat@images$slice1@coordinates[,c("imagerow","imagecol")]
  win <- owin(xrange = c(min(coords$imagecol), max(coords$imagecol)),
              yrange = c(min(coords$imagerow), max(coords$imagerow)))
  Idents(seurat) <- seurat[[col_name]][[col_name]]
  cluster_coord <- coords[WhichCells(seurat, idents = cluster_name),]
  ppp <- ppp(x = cluster_coord$imagecol, y = cluster_coord$imagerow, window = win)
  
  if(!is.null(weights)){
    wts <- seurat@meta.data[WhichCells(seurat, idents = cluster_name),] %>% pull(weights)
    d <- density.ppp(ppp, weights = wts, ...)
  } else {
    d <- density.ppp(ppp,  ...)
  }
  d <- as.function(d)
  density <- d(coords$imagecol, coords$imagerow)
  density <- (density - min(density)) / (max(density) - min(density))
  df <- data.frame(density, row.names = rownames(coords))
  colnames(df) <- paste0(cluster_name, "_density") 
  seurat <- AddMetaData(seurat, df)
  
  return(seurat)
}

generate_cell_density_v2 <- function(seurat, col_name, cluster_name, weights, return_seurat = T, ...){
  library(spatstat)
  library(Seurat)
  
  coords <- GetTissueCoordinates(seurat)
  win <- owin(xrange = c(min(coords$x), max(coords$x)),
              yrange = c(min(coords$y), max(coords$y)))
  Idents(seurat) <- seurat[[col_name]][[col_name]]
  cluster_coord <- coords[WhichCells(seurat, idents = cluster_name),]
  ppp <- ppp(x = cluster_coord$x, y = cluster_coord$y, window = win)
  
  if(!is.null(weights)){
    wts <- seurat@meta.data[WhichCells(seurat, idents = cluster_name),] %>% pull(weights)
    d <- density.ppp(ppp, weights = wts, ...)
  } else {
    d <- density.ppp(ppp,  ...)
  }
  d <- as.function(d)
  density <- d(coords$x, coords$y)
  density <- (density - min(density)) / (max(density) - min(density))
  df <- data.frame(density, row.names = rownames(coords))
  colnames(df) <- paste0(cluster_name, "_density") 
  seurat <- AddMetaData(seurat, df)
  
  if(return_seurat){return(seurat)}else{return(df)}
}

generate_cell_density_v2_multi <- function(seurat, sampleId_col, return_seurat = T, ...){
  pathology_densities <- lapply(unique(seurat[[sampleId_col]][[sampleId_col]]), function(x){
    message('Processing sample: ', x)
    temp_seurat <- seurat[, seurat[[sampleId_col]] == x]
    denisty <- generate_cell_density_v2(temp_seurat, return_seurat = F, ...)
    return(denisty)
  })
  pathology_densities <- bind_rows(pathology_densities)
  seurat <- AddMetaData(seurat, pathology_densities)
  if(return_seurat){return(seurat)}else{return(pathology_densities)}
}


generate_gene_density <- function(seurat, assay, gene, return_seurat = T, ...){
  coords <- GetTissueCoordinates(seurat) 
  win <- owin(xrange = c(min(coords$x), max(coords$x)),
              yrange = c(min(coords$y), max(coords$y)))
  gene_counts <- GetAssayData(seurat, assay = assay, layer = 'data')[gene,]
  pos_cells <- names(gene_counts[gene_counts > 0])
  pos_cells_coord <- coords[pos_cells,]
  ppp <- ppp(x = pos_cells_coord$x, y = pos_cells_coord$y, window = win)
  
  d <- density.ppp(ppp, weights = gene_counts[pos_cells], ...)
  d <- as.function(d)
  density <- d(coords$x, coords$y)
  density <- (density - min(density)) / (max(density) - min(density))
  df <- data.frame(density, row.names = rownames(coords))
  colnames(df) <- paste0(gene, "_density") 
  seurat <- AddMetaData(seurat, df)
  
  if(return_seurat){return(seurat)}else{return(df)}
}

generate_gene_density_multi <- function(seurat, sampleId_col, return_seurat = T, ...){
  gene_densities <- lapply(unique(seurat[[sampleId_col]][[sampleId_col]]), function(x){
    message('Processing sample: ', x)
    temp_seurat <- seurat[, seurat[[sampleId_col]] == x]
    denisty <- generate_gene_density(temp_seurat, return_seurat = F, ...)
    return(denisty)
  })
  gene_densities <- bind_rows(gene_densities)
  seurat <- AddMetaData(seurat, gene_densities)
  if(return_seurat){return(seurat)}else{return(gene_densities)}
}
