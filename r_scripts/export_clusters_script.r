# Get cluster assignments
clusters <- sc@cpart

# Get filtered expression matrix (genes x cells)
fdata <- getfdata(sc)

# Compute cluster centroids (genes x clusters)
cluster_ids <- sort(unique(clusters))
cluster_defs <- sapply(cluster_ids, function(cl) {
  rowMeans(fdata[, clusters == cl, drop = FALSE])
})

# Ensure column names match expected cluster IDs
colnames(cluster_defs) <- as.character(cluster_ids)

# Compute distance of each cell to its cluster centroid
dist_to_centroid <- sapply(names(clusters), function(cell) {
  cl <- as.character(clusters[[cell]])
  sqrt(sum((fdata[, cell] - cluster_defs[, cl])^2))
})

# Cell numbers as sequence
cell_numbers <- seq_along(clusters)

# Combine into dataframe
output <- data.frame(
  CellNumber = cell_numbers,
  Cluster = clusters,
  DistanceToCentroid = dist_to_centroid
)

write.csv(output, "cluster_assignments.csv", row.names = FALSE)

# Save to CSV
write.csv(output, "cluster_assignments.csv", row.names = FALSE)
# Transpose for readability (clusters as rows, genes as columns)
write.csv(t(cluster_defs), "cluster_centroids.csv")

