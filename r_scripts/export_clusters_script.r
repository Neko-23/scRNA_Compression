# Load full gene list from features.tsv (before filtering)
library(Matrix)
library(RaceID)
load("cluster_assignments.RData")
genes <- read.table("../data/GSE294399_WBC_020823_features.tsv/features.tsv", stringsAsFactors = FALSE)[,1]
cells <- read.table("../data/GSE294399_WBC_020823_barcodes.tsv/barcodes.tsv", stringsAsFactors = FALSE)[,1]
# --- Map filtered gene names back to full list ---
all_genes <- genes  # full gene list
all_gene_indices <- setNames(seq_along(all_genes), all_genes)

filtered_genes <- rownames(getfdata(sc))
filtered_gene_indices <- all_gene_indices[filtered_genes]
cell_names <- colnames(getfdata(sc))

# --- 1. Export cluster marker genes using original gene indices ---
all_clusters <- sort(unique(sc@cpart))
marker_genes_per_cluster <- list()

for (cl in all_clusters) {
  de <- clustdiffgenes(sc, cl, pvalue = 0.01)
  upregulated <- rownames(de$dg[de$dg$fc > 1, ])
  marker_indices <- sort(unique(all_gene_indices[upregulated]))
  marker_genes_per_cluster[[as.character(cl)]] <- marker_indices
}

cluster_output <- lapply(names(marker_genes_per_cluster), function(cl) {
  c(as.integer(cl), marker_genes_per_cluster[[cl]])
})

write.table(do.call(rbind, cluster_output),
            file = "../cluster_marker_genes.csv",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

# --- 2. Export cell-cluster gene differences ---
# Build cluster gene sets
cluster_gene_sets <- marker_genes_per_cluster  # already done above

# Create output list
cell_output <- list()

expr_data <- getfdata(sc)

for (cell_idx in seq_along(cell_names)) {
  cell_name <- cell_names[cell_idx]
  cluster <- sc@cpart[cell_name]
  cl_genes <- cluster_gene_sets[[as.character(cluster)]]

  # Get expressed gene indices (in full gene index space)
  expr_values <- expr_data[, cell_name]
  present_filtered <- which(expr_values > 0)
  present_genes <- filtered_gene_indices[present_filtered]

  # Compare to cluster
  present_not_cluster <- setdiff(present_genes, cl_genes)
  absent_but_cluster <- setdiff(cl_genes, present_genes)

  # Format row
  row <- c(cell_idx, cluster, length(present_not_cluster),
           present_not_cluster,
           length(absent_but_cluster),
           absent_but_cluster)
  cell_output[[length(cell_output)+1]] <- row
}

# Pad and write to CSV
max_len <- max(sapply(cell_output, length))
cell_output_padded <- lapply(cell_output, function(row) {
  length(row) <- max_len
  row
})

write.table(do.call(rbind, cell_output_padded),
            file = "../cell_gene_presence.csv",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            na = "")