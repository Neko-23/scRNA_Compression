library(Matrix)
library(RaceID)
load("cluster_assignments.RData")
genes <- read.table("../data/GSE294399_WBC_020823_features.tsv/features.tsv", stringsAsFactors = FALSE)[,1]
cells <- read.table("../data/GSE294399_WBC_020823_barcodes.tsv/barcodes.tsv", stringsAsFactors = FALSE)[,1]

cluster_assignments <- sc@cpart

# Reorder to match the original cell order
ordered_assignments <- cluster_assignments[colnames(getfdata(sc))]

# Write to CSV: each line = cluster number for that cell
write.table(t(ordered_assignments),
            file = "../compressed/raceid/cell_cluster_assignments.csv",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)