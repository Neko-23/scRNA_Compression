library(Matrix)
library(RaceID)
# Read the matrix (replace with the correct path to matrix file)
expr_mat <- readMM("../data/GSE294399_WBC_020823_matrix.mtx/random.mtx")

# Read gene and cell names (replace with the correct paths to gene and cell files)
genes <- read.table("../data/GSE294399_WBC_020823_features.tsv/random.tsv", stringsAsFactors = FALSE)[,1]
cells <- read.table("../data/GSE294399_WBC_020823_barcodes.tsv/random.tsv", stringsAsFactors = FALSE)[,1]

# Assign row and column names
rownames(expr_mat) <- genes
colnames(expr_mat) <- cells

# Convert to SCseq
sc <- SCseq(as.matrix(expr_mat))
# Filter by number of transcripts
sc <- filterdata(sc,mintotal=1)

# Compute distance matrix
sc <- compdist(sc,metric="pearson")

# Perform clustering
sc <- clustexp(sc)

sc <- findoutliers(sc)
sc <- rfcorrect(sc)
sc <- comptsne(sc)
sc <- compfr(sc)
plotmap(sc)
save(sc, file = "cluster_assignments.RData")
