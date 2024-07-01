# Main function to identify rare cells using scSID method
scSID <- function(data, k = 100, h = 0.85, n = 10) {
  
  # Step 1: Identify highly variable genes (Fano genes)
  cat("Step 1: Identifying highly variable genes...\n")
  seurat_obj <- CreateSeuratObject(counts = data)
  seurat_obj <- FindVariableFeatures(object = seurat_obj, selection.method = 'vst', nfeatures = nrow(data), verbose = FALSE)
  vst_values <- seurat_obj@assays$RNA@meta.features$vst.variance.standardized
  density_estimate <- density(vst_values)
  fano_genes <- rownames(data)[vst_values > find_elbow(density_estimate$x[which.max(density_estimate$y):length(density_estimate$x)], 
                                                      density_estimate$y[which.max(density_estimate$y):length(density_estimate$y)])]
  
  # Filter and log transform data for PCA
  filtered_data <- log2(as.matrix(data[fano_genes, ]) + 1)
  
  # Step 2: Perform PCA
  cat("Step 2: Performing PCA...\n")
  pca_result <- irlba(t(filtered_data), nv = min(c(50, nrow(filtered_data) - 1)))
  pca_scores <- t(pca_result$d * t(pca_result$u))
  
  # Step 3: Compute k-nearest neighbors
  cat("Step 3: Computing k-nearest neighbors...\n")
  knn_result <- Neighbour(pca_scores, pca_scores, k = k, build = "kdtree", cores = 0, checks = 1)
  distances <- knn_result$distances
  
  # Step 4: Identify rare cells based on distance gaps
  cat("Step 4: Identifying rare cells based on distance gaps...\n")
  distance_diffs <- distances[, -1, drop = FALSE] - distances[, -ncol(distances), drop = FALSE]
  max_gap_indices <- apply(distance_diffs, 1, which.max)
  
  rare_cells <- list()
  for (l in 2:k) {
    gap_cells <- which(max_gap_indices == l)
    for (cell in gap_cells) {
      neighbors <- knn_result$indices[cell, 1:(l + n)]
      sub_pca <- pca_scores[neighbors, ]
      dist_matrix <- dist(sub_pca)
      hc <- hclust(dist_matrix, method = "single")
      cut_height <- max(hc$height) * h
      cluster_sizes <- table(cutree(hc, h = cut_height))
      
      if (cluster_sizes[1] == l && (hc$height[length(hc$height)] - hc$height[length(hc$height) - 2]) / hc$height[length(hc$height)] > 0.2) {
        rare_cells[[as.character(cell)]] <- knn_result$indices[cell, 1:l]
      }
    }
  }
  
  cat("Rare cell identification complete.\n")
  return(rare_cells)
}
