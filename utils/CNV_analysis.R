## Copied from 
## https://github.com/kharchenkolab/numbat/blob/a367fa55fd3ec6b516c3131b955b61e8d767a722/R/utils.R#L1804
get_clone_profile = function(joint_post, clone_post) {
  joint_post %>%
    inner_join(
      clone_post %>% select(cell, clone = clone_opt),
      by = 'cell'
    ) %>%
    group_by(clone, CHROM, seg, seg_start, seg_end, cnv_state) %>%
    summarise(
      p_cnv = mean(p_cnv),
      size = n(),
      .groups = 'drop'
    ) %>%
    mutate(CHROM = factor(CHROM, 1:22))
}

get_clones_order <- function(metric_mat){
  
  diag(metric_mat) <- 1.0
  metric_mat[lower.tri(metric_mat) & is.na(metric_mat)] <- 0.0
  metric_mat_Dist <- as.dist(1-metric_mat)
  metric_mat[upper.tri(metric_mat)] <- t(metric_mat)[upper.tri(metric_mat)]
  metric_mat[is.na(metric_mat)] <- 0.0
  metric_mat_clust <- hclust(metric_mat_Dist) %>%
    as.dendrogram() %>%
    order.dendrogram()
  order <- rownames(metric_mat)[metric_mat_clust]
  
  return(order)
}
