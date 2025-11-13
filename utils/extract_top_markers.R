top_markers <- function(seurat_markers_table, lfc_cutoff = 0, pct.1_cutoff = 0.2, n = 10){
  seurat_markers_table %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > lfc_cutoff & p_val_adj < 0.05 & pct.1 > pct.1_cutoff) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = n) %>%
    ungroup() %>%
    pull(gene)
}
