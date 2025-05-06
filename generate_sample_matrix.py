import argparse
import numpy as np
import random
import string

parser = argparse.ArgumentParser(
    description="Generate a sample gene-cell expression matrix file."
)
parser.add_argument("--num_cells", type=int, required=True, help="Number of cells")
parser.add_argument("--num_genes", type=int, required=True, help="Number of genes")
parser.add_argument(
    "--output", type=str, default="sample_matrix.txt", help="Output file name"
)
parser.add_argument(
    "--max_expr",
    type=int,
    default=10,
    help="Maximum expression count per gene-cell pair",
)
parser.add_argument("--prob", type=float, default=0.03, help="Probability of each pair")
parser.add_argument("--seed", type=int, default=42, help="Random seed")
args = parser.parse_args()

np.random.seed(args.seed)

# Generate a random sparse matrix (expression counts)
# 3% density of nonzero entries
matrix = np.random.choice(
    [0, 1], size=(args.num_genes, args.num_cells), p=[1-args.prob, args.prob]
) * np.random.randint(1, args.max_expr + 1, size=(args.num_genes, args.num_cells))

with open(args.output, "w") as f:
    f.write("%%MatrixMarket matrix coordinate integer general\n")
    f.write(
        '%metadata_json: {"software_version": "cellranger-7.0.1", "format_version": 2}\n'
    )
    f.write(f"{args.num_genes} {args.num_cells} {(matrix > 0).sum(axis=None)}\n")
    for cell in range(args.num_cells):
        for gene in range(args.num_genes):
            count = matrix[gene, cell]
            if count > 0:
                f.write(f"{gene+1} {cell+1} {count}\n")

# Read barcodes.tsv and features.tsv, and select the first num_cells and num_genes lines
with open(args.output.split("matrix")[0] + "barcodes.tsv/barcodes.tsv", "r") as f:
    barcodes = [line.strip() for _, line in zip(range(args.num_cells), f)]
if len(barcodes) < args.num_cells:
    raise ValueError(f"barcodes.tsv has fewer than {args.num_cells} lines.")
with open(args.output.split("matrix")[0] + "features.tsv/features.tsv", "r") as f:
    features = [line.strip() for _, line in zip(range(args.num_genes), f)]
if len(features) < args.num_genes:
    raise ValueError(f"features.tsv has fewer than {args.num_genes} lines.")
# Overwrite barcodes.tsv and features.tsv with the selected lines
with open(args.output.split("matrix")[0] + "barcodes.tsv/random.tsv", "w") as f:
    for barcode in barcodes:
        f.write(f"{barcode}\n")
with open(args.output.split("matrix")[0] + "features.tsv/random.tsv", "w") as f:
    for feature in features:
        f.write(f"{feature}\n")

print(f"Sample matrix file written to {args.output}")
