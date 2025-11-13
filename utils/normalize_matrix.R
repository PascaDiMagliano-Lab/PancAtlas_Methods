normalize_matrix <- function(mat) {
  min_val <- min(mat)  # Find the minimum value in the matrix
  max_val <- max(mat)  # Find the maximum value in the matrix
  normalized_mat <- (mat - min_val) / (max_val - min_val)  # Apply min-max normalization
  return(normalized_mat)
}
