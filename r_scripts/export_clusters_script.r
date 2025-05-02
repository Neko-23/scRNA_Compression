# Get unique cluster numbers
all_clusters <- sort(unique(sc@cpart))
marker_genes_per_cluster <- list()

# Get gene names and assign indices
gene_names <- rownames(getfdata(sc))
gene_indices <- setNames(seq_along(gene_names), gene_names)

# For each cluster, compute and store sorted, unique marker gene indices
for (cl in all_clusters) {
  de <- clustdiffgenes(sc, cl, pvalue=0.01)
  upregulated <- rownames(de$dg[de$dg$fc > 1, ])
  marker_indices <- sort(unique(gene_indices[upregulated]))
  marker_genes_per_cluster[[as.character(cl)]] <- marker_indices
}

# Format as list of rows for output
output <- lapply(names(marker_genes_per_cluster), function(cl) {
  c(as.integer(cl), marker_genes_per_cluster[[cl]])
})

# Write to CSV
write.table(do.call(rbind, output),
            file = "cluster_marker_genes.csv",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

# part 2: export cell genes

# Get expression data and gene info
expr_mat <- getfdata(sc)
gene_names <- rownames(expr_mat)
gene_indices <- setNames(seq_along(gene_names), gene_names)
cell_names <- colnames(expr_mat)
cell_indices <- setNames(seq_along(cell_names), cell_names)

# Build a gene index list for each cluster (from Part 1)
cluster_gene_sets <- list()
for (cl in all_clusters) {
  de <- clustdiffgenes(sc, cl, pvalue=0.01)
  upregulated <- rownames(de$dg[de$dg$fc > 1, ])
  cluster_gene_sets[[as.character(cl)]] <- gene_indices[upregulated]
}

# Initialize output list
output <- list()

for (cell_name in cell_names) {
  cell_idx <- cell_indices[cell_name]
  cluster <- sc@cpart[cell_name]
  cl_genes <- cluster_gene_sets[[as.character(cluster)]]

  # Get non-zero expression genes for this cell
  expr_values <- expr_mat[, cell_name]
  cell_genes <- which(expr_values > 0)

  # Compare
  present_not_cluster <- setdiff(cell_genes, cl_genes)
  absent_but_cluster <- setdiff(cl_genes, cell_genes)

  # Assemble row
  row <- c(cell_idx, cluster, length(present_not_cluster),
           present_not_cluster,
           length(absent_but_cluster),
           absent_but_cluster)
  output[[length(output)+1]] <- row
}

# Pad rows to equal length
max_len <- max(sapply(output, length))
output_padded <- lapply(output, function(row) {
  length(row) <- max_len
  row
})

# Write to CSV
write.table(do.call(rbind, output_padded), file="cell_gene_presence.csv", sep=",", 
            row.names=FALSE, col.names=FALSE, na="")
