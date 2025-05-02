library(Matrix)
# Read the matrix (replace with the correct path to matrix file)
expr_mat <- readMM("../scRNA_Compression/data/GSE294399_WBC_020823_matrix.mtx/matrix.mtx")

# Read gene and cell names (replace with the correct paths to gene and cell files)
genes <- read.table("../scRNA_Compression/data/GSE294399_WBC_020823_features.tsv/features.tsv", stringsAsFactors = FALSE)[,1]
cells <- read.table("../scRNA_Compression/data/GSE294399_WBC_020823_barcodes.tsv/barcodes.tsv", stringsAsFactors = FALSE)[,1]

# Assign row and column names
rownames(expr_mat) <- genes
colnames(expr_mat) <- cells

library(RaceID)

# Convert to SCseq
sc <- SCseq(as.matrix(expr_mat))
# Filter by number of transcripts
sc <- filterdata(sc,mintotal=2000)

# Compute distance matrix
sc <- compdist(sc,metric="pearson")

# Perform clustering
sc <- clustexp(sc)

sc <- findoutliers(sc)
sc <- rfcorrect(sc)
sc <- comptsne(sc)
sc <- compfr(sc)
plotmap(sc)
