import pickle
from collections import defaultdict
import csv
import numpy as np
from scipy.io import mmread
import heapq
import collections
import os

class Node:
    def __init__(self, symbol=None, frequency=None):
        self.symbol = symbol
        self.frequency = frequency
        self.left = None
        self.right = None

    def __lt__(self, other):
        return self.frequency < other.frequency

def build_huffman_tree(freq_dict):
  
    # Create a priority queue of nodes
    priority_queue = [Node(char, f) for char, f in freq_dict.items()]
    heapq.heapify(priority_queue)

    # Build the Huffman tree
    while len(priority_queue) > 1:
        left_child = heapq.heappop(priority_queue)
        right_child = heapq.heappop(priority_queue)
        merged_node = Node(frequency=left_child.frequency + right_child.frequency)
        merged_node.left = left_child
        merged_node.right = right_child
        heapq.heappush(priority_queue, merged_node)

    return priority_queue[0]

def generate_huffman_codes(node, code="", huffman_codes={}):
    if node is not None:
        if node.symbol is not None:
            huffman_codes[node.symbol] = code
        generate_huffman_codes(node.left, code + "0", huffman_codes)
        generate_huffman_codes(node.right, code + "1", huffman_codes)

    return huffman_codes

def print_tree(node, prefix="", is_left=True):
    if node is not None:
        s = prefix + ("├─0 " if is_left else "└─1 ")
        if node.symbol is not None:
            s += str(node.symbol)
        print(s)
        if node.left or node.right:
            print_tree(node.left, prefix + ("│   " if is_left else "    "), True)
            print_tree(node.right, prefix + ("│   " if is_left else "    "), False)

def compress(cluster_assignments_file, matrix_path, out_path):
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

    matrix = mmread(matrix_path).toarray()
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
    
    if not os.path.isdir(out_path + '/high_level_compress'):
        os.mkdir(out_path + '/high_level_compress')
    with open(out_path + '/high_level_compress/cluster_genes.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        for k in cluster_genes:
            writer.writerow(list(cluster_genes[k]))

    with open(out_path + '/high_level_compress/deltas.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        for d in deltas:
            writer.writerow([cluster_labels[d]] + list(deltas[d]))

    with open(out_path + '/high_level_compress/counts.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        for c in counts:
            writer.writerow(c)

    def low_level_compress(target):
        with open(out_path + '/high_level_compress/' + target + '.csv', newline='') as f:
            lines = f.readlines()

            content = []
            for l in lines:
                content += list(map(int, l.replace('\r', '').replace('\n', '').split(',')))
                content.append('\n')
            content = content[:-1]

            root = build_huffman_tree(collections.Counter(content))
            huffman_codes = generate_huffman_codes(root)

            # Pickle the array and save to a file
            if not os.path.isdir(out_path + '/low_level_compress'):
                os.mkdir(out_path + '/low_level_compress')
            filename = out_path + '/low_level_compress/' + target + '_huffman_tree'
            with open(filename, 'wb') as file:
                pickle.dump((root, len(content)), file)
            file.close()
            
            with open(out_path + '/low_level_compress/huffman_encoded_'+target, 'wb') as file:
                leftover = ''
                for index, l in enumerate(lines):
                    elems = l.replace('\r', '').replace('\n', '').split(',')
                    elems = list(map(int, elems))
                    if index < len(lines) - 1:
                        elems.append('\n')
                    binary_string = leftover
                    for e in elems:
                        binary_string += huffman_codes[e]

                    main_part = binary_string[:(8 * (len(binary_string) // 8))]
                    leftover = binary_string[(8 * (len(binary_string) // 8)):]

                    for i in range(0, len(main_part), 8):
                        byte_string = main_part[i:i + 8]
                        byte_value = int(byte_string, 2)
                        file.write(bytes([byte_value]))

                if leftover:
                    leftover = leftover.ljust(8, '0')
                    byte_value = int(leftover, 2)
                    file.write(bytes([byte_value]))
            file.close()
        f.close()

    low_level_compress('deltas')
    low_level_compress('counts')

def high_level_decompress(in_path, matrix_path):
    gene_map = defaultdict(lambda: [None, set()])
    cluster_map = defaultdict(lambda: set())
    counts = []
    with open(in_path + '/high_level_compress/deltas.csv', newline='') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            cluster = int(row[0])
            genes = set(map(int, row[1:])) if len(row) > 1 else set()
            gene_map[i] = [cluster, genes]
            cluster_map[cluster].add(i)
            num_cells = i

    with open(in_path + '/high_level_compress/cluster_genes.csv', newline='') as f:
        reader = csv.reader(f)
        for cluster, row in enumerate(reader):
            genes = set(map(int, row))
            for cell in cluster_map[cluster]:
                gene_map[cell][1] = gene_map[cell][1].union(genes)

    with open(in_path + '/high_level_compress/counts.csv', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            counts.append(row)

    matrix = mmread(matrix_path).toarray()
    cell_gene_matrix = np.zeros(matrix.shape)
    for cell_idx in range(0, num_cells+1):
        for i, g in enumerate(sorted(gene_map[cell_idx][1])):
            cell_gene_matrix[g, cell_idx] = counts[cell_idx][i]
    
    print("Accuracy Check " + "Passed" if (cell_gene_matrix != matrix).sum(axis=None) == 0 else "Failed")