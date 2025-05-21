# Single-cell RNA Compression using Clusters and Huffman Encoding

## Abstract
Single-cell RNA sequencing (scRNA-seq) technologies produce large and sparse count matrices that pose significant storage challenges, especially when scaling to datasets containing many samples. This project aims to develop a lossless delta encoding method to substantially reduce the storage footprint of scRNA count matrices compressed using current formats such as CSR, CSC, MTX, and gzip. Our approach leverages existing clustering algorithms to identify groups of similar cells based on their gene expression profiles. For each cluster, we generate a set of reference genes composed of genes commonly expressed across all cells in the cluster. Each individual cell is then represented by storing only genes that differ from the reference in the cluster. The delta-compressed files are then further compressed using Huffman encoding. We evaluated the effectiveness of our method by applying it to a dataset comprising multiple scRNA count matrices and compared the storage of our compressed file with other compression formats. We found that our compression method reduced the storage space by over 83 percent for one of the scRNA matrices stored MTX format.

## Running our Compression Pipeline
1. Retrieve cluster labels by running k-means or our neural.ipynb notebook.
2. Set the sample dataset number and path that you wish to store the output in compression_pipeline.ipynb
3. Run the compression_pipeline.ipynb notebook
# Single-cell RNA Compression using Clusters and Huffman Encoding

## Overview

Single-cell RNA sequencing (scRNA-seq) technologies produce large and sparse count matrices that pose significant storage challenges, especially when scaling to datasets containing many samples. This project develops a lossless delta encoding method to substantially reduce the storage footprint of scRNA count matrices compressed using current formats such as CSR, CSC, MTX, and gzip.

Our approach leverages clustering algorithms to identify groups of similar cells based on their gene expression profiles. For each cluster, we generate a set of reference genes composed of genes commonly expressed across all cells in the cluster. Each individual cell is then represented by storing only genes that differ from the reference in the cluster. The delta-compressed files are then further compressed using Huffman encoding.

We evaluated the effectiveness of our method by applying it to a dataset comprising multiple scRNA count matrices and compared the storage of our compressed file with other compression formats. We found that our compression method reduced the storage space by over 83 percent for one of the scRNA matrices stored in MTX format.

All code is available at the [GitHub repository](https://github.com/Neko-23/scRNA_Compression).

---

## Table of Contents

- [Introduction](#introduction)
- [Related Work](#related-work)
- [Methods](#methods)
  - [Data Acquisition and Preprocessing](#data-acquisition-and-preprocessing)
  - [Cluster-Based Compression Pipeline](#cluster-based-compression-pipeline)
  - [High-level Compression with Clusters and Delta Encoding](#high-level-compression-with-clusters-and-delta-encoding)
  - [Low-level Compression with Huffman Encoding](#low-level-compression-with-huffman-encoding)
  - [Evaluation and Benchmarking](#evaluation-and-benchmarking)
  - [Clustering Methods](#clustering-methods)
- [Results](#results)
- [Conclusion and Discussion](#conclusion-and-discussion)
- [References](#references)

---

## Introduction

Single-cell RNA sequencing (scRNA-seq) has revolutionized genomics by enabling researchers to measure gene expression at the single-cell level. Unlike bulk RNA sequencing, which provides an averaged expression profile, scRNA-seq captures the heterogeneity of gene expression between individual cells. This makes scRNA-seq a powerful tool for studying complex biological systems, such as cellular differentiation, tissue composition, and disease progression.

The resulting gene expression matrices are extremely sparse, with most entries being zero. While formats like MTX, CSR, and CSC efficiently store sparse matrices, further compression is possible by leveraging biological structure, such as cell clustering.

---

## Related Work

- **Sparse Matrix Formats:** MTX, CSR, and CSC are standard formats for storing sparse matrices. They reduce storage by only recording nonzero entries and compressing row/column indices.
- **Clustering in scRNA-seq:** Clustering is widely used for cell taxonomy and function identification. Libraries such as SAIC, RaceID, and scDeepCluster provide various clustering algorithms tailored for scRNA-seq data.
- **Lossless Compression:** Huffman encoding and related schemes (e.g., Burrows-Wheeler transform, move-to-front, run-length encoding) are established lossless compression techniques. Recent work like ScBlkCom applies such methods specifically to scRNA-seq data.

---

## Methods

### Data Acquisition and Preprocessing

- Two representative scRNA-seq datasets (Sample 1 and Sample 2) were obtained from the Gene Expression Omnibus.
- Preprocessing included quality control, normalization, and conversion to MTX, CSR, and CSC formats using Python (`scipy`, `numpy`).
- Both compressed and uncompressed versions of each format were benchmarked.

### Cluster-Based Compression Pipeline

- Cells are clustered using k-means or a neural clustering approach.
- Cluster labels are used for high-level compression (delta encoding), followed by low-level compression (Huffman encoding).
- All code and scripts are available in this repository.

### High-level Compression with Clusters and Delta Encoding

- For each cluster, a reference gene set is constructed from genes commonly expressed in all cluster cells.
- Each cell is represented by the genes that differ from the cluster reference (delta).
- Data is stored in three files: `cluster_genes.csv`, `deltas.csv`, and `counts.csv`.

### Low-level Compression with Huffman Encoding

- Huffman encoding is applied to `deltas.csv` and `counts.csv` (not to `cluster_genes.csv` due to its small size).
- Huffman trees are stored as pickled `.pkl` objects.
- Both uncompressed and gzip-compressed Huffman-encoded files are produced.

### Evaluation and Benchmarking

- File sizes are measured using Unix `du` and Python file I/O.
- Compression ratios are calculated as the ratio of compressed to original MTX file size.
- Only file sizes are measured; runtime and memory usage are not the focus.

### Clustering Methods

- **K-means:** Standard k-means from `scikit-learn` (random seed 42), with varying $k$.
- **Neural Clustering:** Autoencoder-based clustering, followed by k-means in latent space.
- **RaceID:** Specialized R library for rare cell type identification (used for comparison).

---

## Results

- **Standard Formats:** MTX, CSR, and CSC (with and without gzip) yield file sizes between 7.75 MB and 57.12 MB.
- **Cluster-Based Compression:** Both k-means and neural clustering methods yield substantial reductions in file size, especially with delta and Huffman encoding.
- **Best Results:** Neural clustering with delta and Huffman encoding plus gzip achieves the smallest file sizes (e.g., 7.64 MB for Sample 1).
- **Cluster Statistics:** Detailed statistics for each clustering method and $k$ value are provided in the manuscript and supplemental tables.
- **Compression Ratios:** The best methods achieve over 83% reduction in file size compared to the original MTX format.

---

## Conclusion and Discussion

Our cluster-based compression pipeline for scRNA-seq count matrices, combining clustering, delta encoding, and Huffman encoding, yields substantial storage reductions compared to standard formats. Both k-means and neural clustering approaches provide consistent gains, with neural clustering offering additional benefits for more complex datasets.

The approach is lossless and preserves the full information content, ensuring compatibility with downstream analyses. While clustering and reference construction introduce computational overhead, the resulting storage savings are significant, especially as single-cell datasets grow in scale.

**Future Work:**

- Explore adaptive clustering strategies and alternative neural architectures.
- Benchmark computational overhead (runtime, memory).
- Extend the pipeline to other types of high-dimensional, sparse biological data.
- Incorporate additional compression steps (e.g., Burrows-Wheeler transform, run-length encoding).

---

## References

Please see the `reference.bib` file for all citations, including:

- Prior work on scRNA-seq clustering and compression
- Libraries and tools used (scikit-learn, RaceID, GeeksforGeeks Huffman code, etc.)
- Datasets and algorithms referenced in the manuscript

---

## Source Code

All code for this project is available at:  
[https://github.com/Neko-23/scRNA_Compression](https://github.com/Neko-23/scRNA_Compression)

If the repository is private, please ensure that all relevant reviewers have access.

---

## Citation

If you use this code or method, please cite our manuscript and the relevant prior work as described in the references.

---
