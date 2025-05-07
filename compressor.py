import pickle
from collections import defaultdict
import csv
import numpy as np

def compress(cluster_assignments_file, cluster_out_file, deltas_file, matrix):
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
    for i in range(matrix.shape[1]):# Cell
        genes_set = set(np.arange(0, matrix.shape[0] + 0)[matrix[:,i].astype(bool)].tolist())
        deltas[i] = genes_set - cluster_genes[cluster_labels[i]]
        
    with open(cluster_out_file, 'w', newline='') as file:
        writer = csv.writer(file)
        for k in cluster_genes:
            writer.writerow([k] + list(cluster_genes[k]))

    with open(deltas_file, 'w', newline='') as file:
        writer = csv.writer(file)
        for d in deltas:
            writer.writerow([cluster_labels[d]] + list(deltas[d]))

def decompress(clusters_file, deltas_file, shape):
    gene_map = defaultdict(lambda: [None, set()])
    cluster_map = defaultdict(lambda: set())
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
        for row in reader:
            cluster = int(row[0])
            genes = set(map(int, row[1:]))
            for cell in cluster_map[cluster]:
                gene_map[cell][1] = gene_map[cell][1].union(genes)

    cell_gene_matrix = np.zeros(shape)
    for cell_idx in range(0, num_cells+1):
        for g in gene_map[cell_idx][1]:
            cell_gene_matrix[g, cell_idx] = 1
    
    return cell_gene_matrix