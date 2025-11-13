detect_singlets <- function(seurat, col, ...){
  coords <- GetTissueCoordinates(seurat)
  win <- owin(xrange = c(min(coords$x), max(coords$x)),
              yrange = c(min(coords$y), max(coords$y)))
  Idents(seurat) <- seurat[[col]][[col]]
  
  clusters <- unique(as.character(seurat@meta.data[[col]]))
  cells_class <- lapply(clusters, function(x){
    cluster_coord <- coords[WhichCells(seurat, idents = x),]
    ppp <- ppp(x = cluster_coord$x, y = cluster_coord$y, window = win, marks = cluster_coord)
    ppp <- nnclean(ppp, ...)
    ppp$marks %>% select(class, prob)
  }) %>% 
    bind_rows()
  
  return(cells_class)
}

detect_boundary <- function(seurat, col){
  seurat_coord <- GetTissueCoordinates(seurat)
  win <- owin(xrange = c(min(seurat_coord$x), max(seurat_coord$x)),
              yrange = c(min(seurat_coord$y), max(seurat_coord$y)))
  ann <- seurat@meta.data[[col]]
  ppp <- ppp(x = seurat_coord$x, y = seurat_coord$y,
             window = win,
             marks = data.frame(cellname = colnames(seurat),
                                annotation = seurat@meta.data[[col]]))
  
  nn <- nnwhich(ppp, k = c(1,2,3,4,5,6)) %>%
    data.frame() %>%
    `rownames<-`(colnames(seurat))
  
  nn <- apply(nn, 1, function(x){return(ppp$marks$annotation[x])}) %>%
    t() %>%
    data.frame() %>%
    rownames_to_column(var = 'index_cell')
  
  nn$index_cellType <- merge(nn %>% select(index_cell),
                             ppp$marks, 
                             by.x = 'index_cell',
                             by.y = 'cellname') %>% 
    pull(annotation)
  
  nn <- column_to_rownames(nn, 'index_cell')
  
  nhood_type <- apply(nn, 1, function(x){
    ifelse(length(unique(x)) > 1, "Heterogenic", "Homogenic")
  })
  
  nn[[paste0(col, "_nhood_type")]] <- nhood_type
  # nn$cell_annotation <- seurat@meta.data[[col]]
  # nn$nhood_CT <- paste0(nn[[paste0(col, "_nhood_type")]], "_", nn$cell_annotation)
  seurat <- AddMetaData(seurat, nn %>% select(paste0(col, "_nhood_type")))
  return(seurat)
}

calculate_JSD <- function(seurat, col, cluster1, cluster2, ...){
  library(philentropy)
  # extracting coordinates
  coords <- GetTissueCoordinates(seurat)
  win <- owin(xrange = c(min(coords$x), max(coords$x)),
              yrange = c(min(coords$y), max(coords$y)))
  Idents(seurat) <- seurat[[col]][[col]]
  
  # density for first cluster
  cluster_coord <- coords[WhichCells(seurat, idents = cluster1),]
  ppp <- ppp(x = cluster_coord$x, y = cluster_coord$y, window = win)
  d1 <- density.ppp(ppp, ...)
  d1_vec <- c(d1$v)
  d1_vec <- (d1_vec - min(d1_vec)) / (max(d1_vec) - min(d1_vec))
  
  # density for second cluster
  cluster_coord <- coords[WhichCells(seurat, idents = cluster2),]
  ppp <- ppp(x = cluster_coord$x, y = cluster_coord$y, window = win)
  d2 <- density.ppp(ppp, sigma = bw.diggle)
  d2_vec <- c(d2$v)
  d2_vec <- (d2_vec - min(d2_vec)) / (max(d2_vec) - min(d2_vec))
  
  # calculate JSD
  jsd <- JSD(rbind(d1_vec, d2_vec), est.prob = 'empirical')
  
  return(jsd)
}

pairwise_JSD <- function(seurat, col, neg_log = TRUE){
  clusters <- as.character(unique(seurat@meta.data[[col]]))
  JSD_pairwise_distance <- lapply(clusters, function(x){
    distances <- lapply(clusters, function(y){
      jsd <- calculate_JSD(seurat, col, x, y, sigma = bw.diggle)
    }) %>% 
      `names<-`(clusters) %>%
      unlist() %>%
      data.frame() %>%
      rownames_to_column() %>%
      mutate(cluster1 = x,
             rowname = gsub('.jensen-shannon','', rowname)) %>%
      `colnames<-`(c('cluster2', 'JSD', 'cluster1')) %>%
      select(cluster1, cluster2, JSD)
  }) %>% bind_rows()
  
  if(neg_log){
  JSD_pairwise_distance$JSD <- -log(JSD_pairwise_distance$JSD) 
  JSD_pairwise_distance$JSD[is.infinite(JSD_pairwise_distance$JSD)] <- NA
  }
  
  JSD_pairwise_distance <- pivot_wider(JSD_pairwise_distance, 
                                       names_from = cluster2, 
                                       values_from = JSD) %>%
    column_to_rownames('cluster1')
    return(JSD_pairwise_distance)
}

compute_average_min_distance <- function(seurat, col, log = T){
  coords <- GetTissueCoordinates(seurat)
  win <- owin(xrange = c(min(coords$x), max(coords$x)),
              yrange = c(min(coords$y), max(coords$y)))
  Idents(seurat) <- seurat[[col]][[col]]
  
  clusters <- as.character(unique(seurat@meta.data[[col]]))

  pairwise_distance <- lapply(clusters, function(x){
    distances <- lapply(clusters, function(y){
      cluster1_coord <- coords[WhichCells(seurat, idents = x),]
      cluster2_coord <- coords[WhichCells(seurat, idents = y),]
      
      cluster1_ppp <- ppp(x = cluster1_coord$x, y = cluster1_coord$y, window = win)
      cluster2_ppp <- ppp(x = cluster2_coord$x, y = cluster2_coord$y, window = win)
      
      dist <- crossdist(cluster1_ppp, cluster2_ppp) 
      return(mean(rowMins(dist)))
      # return(mean(dist))
    }) %>% `names<-`(clusters) %>% unlist() %>%
      data.frame() %>% rownames_to_column() %>%
      mutate(cluster1 = x) %>% `colnames<-`(c("cluster2", "avg_distance", "cluster1")) %>%
      select(cluster1, cluster2, avg_distance)
  }) %>% bind_rows()
  
  if(log){
    pairwise_distance$avg_distance <- log(pairwise_distance$avg_distance)
    pairwise_distance$avg_distance[is.infinite(pairwise_distance$avg_distance)] <- NA
    # pairwise_distance$avg_distance <- scale(pairwise_distance$avg_distance)
  }
  
  pairwise_distance <- pivot_wider(pairwise_distance, 
                                   names_from = cluster2, 
                                   values_from = avg_distance) %>%
    column_to_rownames('cluster1')
  
  pairwise_distance[is.na(pairwise_distance)] <- 0
  
  return(pairwise_distance)
}
