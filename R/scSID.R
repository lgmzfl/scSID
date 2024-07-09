# Function to identify rare cells using scSID methodology
scSID <- function(data, k = 100, h = 0.85, n = 10) {
  
  # Step 1: Identification of Fano genes for clustering purposes
  cat("Step 1: Identifying highly variable (Fano) genes...\n")
  seurat_obj <- CreateSeuratObject(counts = data)
  obj <- FindVariableFeatures(object = seurat_obj, selection.method = 'vst', nfeatures = nrow(data), verbose = FALSE)
  vst_values <- obj@assays$RNA@meta.features$vst.variance.standardized
  density_estimation <- density(vst_values)
  fano_genes <- rownames(data)[vst_values > elbow(density_estimation$x[which.max(density_estimation$y):length(density_estimation$x)], 
                                                            density_estimation$y[which.max(density_estimation$y):length(density_estimation$y)])]
  
  # Log transformation of the filtered data
  filtered_data <- log2(as.matrix(data[fano_genes, ]) + 1)
  
  # Step 2: Principal Component Analysis (PCA)
  cat("Step 2: Performing PCA...\n")
  pca_result <- irlba(t(filtered_data), nv = min(c(50, nrow(filtered_data) - 1)))
  pca_scores <- t(pca_result$d * t(pca_result$u))
  
  # Step 3: K-Nearest Neighbors computation
  cat("Step 3: Computing K-Nearest Neighbors...\n")
  knn_result <- Neighb(pca_scores, pca_scores, k = k, build = "kdtree", cores = 0, checks = 1)
  distances <- knn_result$distances
  
  # Step 4: Identification of rare cells based on distance differentials
  cat("Step 4: Identifying rare cells based on distance differentials...\n")
  distance_differentials <- distances[, -1, drop = FALSE] - distances[, -ncol(distances), drop = FALSE]
  max_gap_indices <- apply(distance_differentials, 1, which.max)
  
  # Call the identify_rare_cells function
  rare_cells <- identify_rare_cells(k, max_gap_indices, knn_result, pca_scores, h, n)
  
  cat("Rare cell identification complete.\n")
  return(rare_cells)
}