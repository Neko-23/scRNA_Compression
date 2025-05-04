# Load full gene list from features.tsv (before filtering)
all_genes <- read.table("../data/GSE294399_WBC_020823_features.tsv/features.tsv", stringsAsFactors = FALSE)[,1]
all_gene_indices <- setNames(seq_along(all_genes), all_genes)

# Get all unique cluster IDs
all_clusters <- sort(unique(sc@cpart))
marker_genes_per_cluster <- list()

# Compute marker genes using full indexing
for (cl in all_clusters) {
  de <- clustdiffgenes(sc, cl, pvalue = 0.01)
  upregulated <- rownames(de$dg[de$dg$fc > 1, ])
  marker_indices <- sort(unique(all_gene_indices[upregulated]))
  marker_genes_per_cluster[[as.character(cl)]] <- marker_indices
}

# Format for CSV
output <- lapply(names(marker_genes_per_cluster), function(cl) {
  c(as.integer(cl), marker_genes_per_cluster[[cl]])
})

# Save to CSV
write.table(do.call(rbind, output),
            file = "../cluster_marker_genes.csv",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

# Load full gene list again
all_genes <- read.table("../data/GSE294399_WBC_020823_features.tsv/features.tsv", stringsAsFactors = FALSE)[,1]
all_gene_indices <- setNames(seq_along(all_genes), all_genes)

# Prepare cell and gene data
expr_mat <- getfdata(sc)  # filtered matrix
gene_names_filtered <- rownames(expr_mat)
cell_names <- colnames(expr_mat)

# Map filtered gene names to full indices
filtered_gene_indices <- all_gene_indices[gene_names_filtered]

# Build cluster gene sets using full indexing
all_clusters <- sort(unique(sc@cpart))
cluster_gene_sets <- list()

for (cl in all_clusters) {
  de <- clustdiffgenes(sc, cl, pvalue = 0.01)
  upregulated <- rownames(de$dg[de$dg$fc > 1, ])
  cluster_gene_sets[[as.character(cl)]] <- sort(unique(all_gene_indices[upregulated]))
}

# Generate output rows
output <- list()

for (cell_idx in seq_along(cell_names)) {
  cell_name <- cell_names[cell_idx]
  cluster <- sc@cpart[cell_name]
  cl_genes <- cluster_gene_sets[[as.character(cluster)]]

  # Get expressed gene indices (based on original indices)
  expr_values <- expr_mat[, cell_name]
  present_filtered <- which(expr_values > 0)
  present_genes <- filtered_gene_indices[present_filtered]

  # Compare sets
  present_not_cluster <- setdiff(present_genes, cl_genes)
  absent_but_cluster <- setdiff(cl_genes, present_genes)

  # Assemble output row
  row <- c(cell_idx, cluster, length(present_not_cluster),
           present_not_cluster,
           length(absent_but_cluster),
           absent_but_cluster)
  output[[length(output)+1]] <- row
}

# Pad rows for CSV
max_len <- max(sapply(output, length))
output_padded <- lapply(output, function(row) {
  length(row) <- max_len
  row
})

# Write to CSV
write.table(do.call(rbind, output_padded),
            file = "../cell_gene_presence.csv",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            na = "")