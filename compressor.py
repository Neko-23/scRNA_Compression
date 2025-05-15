import pickle
from collections import defaultdict
import csv
import numpy as np

def compress(cluster_assignments_file, cluster_out_file, deltas_file, counts_file, matrix):
    # Load the pickled array from the file
    extension = cluster_assignments_file.split('.')[-1]
    if extension == 'pkl':
        with open(cluster_assignments_file, 'rb') as file:
            cluster_labels = pickle.load(file)
    elif extension == 'csv':
        cluster_labels = []
        with open(cluster_assignments_file, newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                if row:
                    cluster_labels += list(map(lambda a: int(a) if a else None,row))
    else:
        raise NotImplementedError(f"Unsupported file type: *.{extension}")

    genes = defaultdict(lambda: [])
    cluster_genes = defaultdict(lambda: set())
    for i in range(matrix.shape[1]):# Cell
        genes[cluster_labels[i]].append(set(np.arange(0, matrix.shape[0] + 0)[matrix[:,i].astype(bool)].tolist()))
    for i in range(max(cluster_labels) + 1):
        if genes[i]:
            cluster_genes[i] = set.intersection(*genes[i])
            
    deltas = defaultdict(lambda: set())
    counts = []
    for i in range(matrix.shape[1]):# Cell
        genes_set = set(np.arange(0, matrix.shape[0] + 0)[matrix[:,i].astype(bool)].tolist())
        deltas[i] = genes_set - cluster_genes[cluster_labels[i]]
        counts.append(matrix[:,i][matrix[:,i] > 0].tolist())
        
    with open(cluster_out_file, 'w', newline='') as file:
        writer = csv.writer(file)
        for k in cluster_genes:
            writer.writerow(list(cluster_genes[k]))

    with open(deltas_file, 'w', newline='') as file:
        writer = csv.writer(file)
        for d in deltas:
            writer.writerow([cluster_labels[d]] + list(deltas[d]))

    with open(counts_file, 'w', newline='') as file:
        writer = csv.writer(file)
        for c in counts:
            writer.writerow(c)

def decompress(clusters_file, deltas_file, counts_file, shape):
    gene_map = defaultdict(lambda: [None, set()])
    cluster_map = defaultdict(lambda: set())
    counts = []
    with open(deltas_file, newline='') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            cluster = int(row[0])
            genes = set(map(int, row[1:])) if len(row) > 1 else set()
            gene_map[i] = [cluster, genes]
            cluster_map[cluster].add(i)
            num_cells = i

    with open(clusters_file, newline='') as f:
        reader = csv.reader(f)
        for cluster, row in enumerate(reader):
            genes = set(map(int, row))
            for cell in cluster_map[cluster]:
                gene_map[cell][1] = gene_map[cell][1].union(genes)

    with open(counts_file, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            counts.append(row)

    cell_gene_matrix = np.zeros(shape)
    for cell_idx in range(0, num_cells+1):
        for i, g in enumerate(sorted(gene_map[cell_idx][1])):
            cell_gene_matrix[g, cell_idx] = counts[cell_idx][i]
    
    return cell_gene_matrix