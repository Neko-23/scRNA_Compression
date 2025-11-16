import pickle
from collections import defaultdict
import csv
import numpy as np
from scipy.io import mmread
from scipy.sparse import coo_matrix, csr_matrix
import heapq
import collections
import os
import shutil

#fpgrowth libraries
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth
import pandas as pd

# Huffman code retrieved from https://www.geeksforgeeks.org/huffman-coding-in-python/
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

# helper copy file function
def copy_file(source_path, destination_path):
    """Copies a file from the source path to the destination path.

    Args:
        source_path: The path to the file to be copied.
        destination_path: The path to the destination directory.
    """
    try:
        shutil.copy(source_path, destination_path)
        print(f"File copied successfully from {source_path} to {destination_path}")
    except FileNotFoundError:
        print(f"Error: Source file not found: {source_path}")
    except PermissionError:
        print(f"Error: Permission denied to access {source_path} or {destination_path}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

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

    high_level_dir = os.path.join(out_path, "high_level_compress")
    def low_level_compress(target):
        src_file = os.path.join(high_level_dir, f"{target}.csv")
        with open(src_file, newline="") as f:
            lines = f.readlines()

        content = []
        for l in lines:
            # tolerate stray whitespace and empty tokens
            parts = l.replace("\r", "").replace("\n", "").split(",")
            for p in parts:
                p = p.strip()
                if p == "":
                    continue
                try:
                    content.append(int(p))
                except ValueError:
                    # ignore malformed tokens
                    continue
            content.append("\n")
        if content:
            content = content[:-1]

        root = build_huffman_tree(collections.Counter(content))
        huffman_codes = generate_huffman_codes(root)

        # Pickle the array and save to a file
        low_level_dir = os.path.join(out_path, "low_level_compress")
        os.makedirs(low_level_dir, exist_ok=True)
        filename = os.path.join(low_level_dir, f"{target}_huffman_tree")
        with open(filename, "wb") as file:
            pickle.dump((root, len(content)), file)

        with open(
            os.path.join(low_level_dir, f"huffman_encoded_{target}"), "wb"
        ) as file:
            leftover = ""
            for index, l in enumerate(lines):
                elems = [
                    p
                    for p in l.replace("\r", "").replace("\n", "").split(",")
                    if p.strip() != ""
                ]
                try:
                    elems = list(map(int, elems))
                except ValueError:
                    # skip malformed lines
                    elems = []
                if index < len(lines) - 1:
                    elems.append("\n")
                binary_string = leftover
                for e in elems:
                    binary_string += huffman_codes[e]

                main_part = binary_string[: (8 * (len(binary_string) // 8))]
                leftover = binary_string[(8 * (len(binary_string) // 8)) :]

                for i in range(0, len(main_part), 8):
                    byte_string = main_part[i : i + 8]
                    byte_value = int(byte_string, 2)
                    file.write(bytes([byte_value]))

            if leftover:
                leftover = leftover.ljust(8, "0")
                byte_value = int(leftover, 2)
                file.write(bytes([byte_value]))

    low_level_compress("deltas")
    low_level_compress("counts")

def high_level_decompress(in_path, matrix_path):
    gene_map = defaultdict(lambda: [None, set()])
    cluster_map = defaultdict(lambda: set())
    counts = []
    with open(in_path + '/deltas.csv', newline='') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            cluster = int(row[0])
            genes = set(map(int, row[1:])) if len(row) > 1 else set()
            gene_map[i] = [cluster, genes]
            cluster_map[cluster].add(i)
            num_cells = i
        f.close()

    with open(in_path + '/cluster_genes.csv', newline='') as f:
        reader = csv.reader(f)
        for cluster, row in enumerate(reader):
            genes = set(map(int, row))
            for cell in cluster_map[cluster]:
                gene_map[cell][1] = gene_map[cell][1].union(genes)

    with open(in_path + '/counts.csv', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            counts.append(row)

    matrix = mmread(matrix_path).toarray()
    cell_gene_matrix = np.zeros(matrix.shape)
    for cell_idx in range(0, num_cells+1):
        for i, g in enumerate(sorted(gene_map[cell_idx][1])):
            cell_gene_matrix[g, cell_idx] = counts[cell_idx][i]
    
    print("High Level Accuracy Check " + "Passed" if (cell_gene_matrix != matrix).sum(axis=None) == 0 else "Failed")    

def low_level_decompress(in_path):
    def helper_decompress(target):
        tree_file = os.path.join(
            in_path, "low_level_compress", f"{target}_huffman_tree"
        )
        with open(tree_file, "rb") as file:
            # Load the pickled object from the file
            root, length = pickle.load(file)

        decoded_lines = []
        count = 0
        enc_file = os.path.join(
            in_path, "low_level_compress", f"huffman_encoded_{target}"
        )
        with open(enc_file, "rb") as file:
            decoded_line = []
            cur = root
            stop = False
            while byte := file.read(1):
                b = byte[0]
                for i in range(7, -1, -1):
                    if count >= length:
                        stop = True
                        break
                    # traverse the tree according to the next bit
                    bit = (b >> i) & 1
                    cur = cur.right if bit else cur.left
                    if cur is None:
                        # malformed traversal or padding: reset
                        cur = root
                        continue
                    # when we reach a leaf, emit symbol
                    if cur.symbol is not None:
                        count += 1
                        sym = cur.symbol
                        if sym != "\n":
                            decoded_line.append(sym)
                        else:
                            decoded_lines.append(decoded_line)
                            decoded_line = []
                        cur = root
                if stop:
                    break
            if decoded_line:
                decoded_lines.append(decoded_line)

        result = "ENCODED CORRECTLY!"
        hl_file = os.path.join(in_path, "high_level_compress", f"{target}.csv")
        with open(hl_file, newline="") as f:
            lines = f.readlines()
        for i, l in enumerate(lines):
            parts = [
                p.strip()
                for p in l.replace("\r", "").replace("\n", "").split(",")
                if p.strip() != ""
            ]
            try:
                row_ints = list(map(int, parts))
            except ValueError:
                row_ints = []
            if i >= len(decoded_lines) or row_ints != decoded_lines[i]:
                result = "ENCODED WRONG!!"
                break
        print(target + " huffman : " + result)

    helper_decompress("deltas")
    helper_decompress("counts")