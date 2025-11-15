import pickle
from collections import defaultdict
import csv
import numpy as np
from scipy.io import mmread
from scipy.sparse import coo_matrix, csr_matrix
import heapq
import collections
import os


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


def compress(cluster_assignments_file, matrix_path, out_path):
    # Load the pickled array from the file
    extension = cluster_assignments_file.split(".")[-1]
    if extension == "pkl":
        with open(cluster_assignments_file, "rb") as file:
            cluster_labels = pickle.load(file)
    elif extension == "csv":
        cluster_labels = []
        with open(cluster_assignments_file, newline="") as f:
            reader = csv.reader(f)
            for row in reader:
                if row:
                    cluster_labels += list(map(lambda a: int(a) if a else None, row))
    else:
        raise NotImplementedError(f"Unsupported file type: *.{extension}")

    # operate on sparse matrix to avoid dense conversions
    sparse_matrix = mmread(matrix_path)
    # use CSC for efficient column access
    sparse_csc = sparse_matrix.tocsc()
    G, N = sparse_csc.shape
    genes = defaultdict(lambda: [])
    cluster_genes = defaultdict(lambda: set())
    for i in range(N):  # Cell
        start = sparse_csc.indptr[i]
        end = sparse_csc.indptr[i + 1]
        rows = sparse_csc.indices[start:end]
        genes[cluster_labels[i]].append(set(rows.tolist()))
    for i in range(max(cluster_labels) + 1):
        if genes[i]:
            cluster_genes[i] = set.intersection(*genes[i])

    deltas = defaultdict(lambda: set())
    counts = []
    compressed_counts = []
    for i in range(N):  # Cell
        start = sparse_csc.indptr[i]
        end = sparse_csc.indptr[i + 1]
        rows = sparse_csc.indices[start:end]
        vals = sparse_csc.data[start:end]
        genes_set = set(rows.tolist())
        deltas[i] = genes_set - cluster_genes[cluster_labels[i]]
        # ensure counts are ordered by sorted gene indices for reproducibility
        if len(rows) > 0:
            order = np.argsort(rows)
            sorted_vals = [int(v) for v in vals[order]]
            sorted_rows = rows[order]
        else:
            sorted_vals = []
            sorted_rows = []
        next_row = sorted_vals
        next_row_compressed = []
        counter = 0
        # use a different loop variable name and emit trailing run token if present
        for val in next_row:
            if val == 1:
                counter += 1
            else:
                if counter:
                    next_row_compressed.append(-counter)
                    counter = 0
                next_row_compressed.append(int(val))
        if counter:
            next_row_compressed.append(-counter)
        counts.append(next_row)
        compressed_counts.append(next_row_compressed)

    high_level_dir = os.path.join(out_path, "high_level_compress_count_extension")
    os.makedirs(high_level_dir, exist_ok=True)
    with open(
        os.path.join(high_level_dir, "cluster_genes.csv"), "w", newline=""
    ) as file:
        writer = csv.writer(file)
        for k in cluster_genes:
            writer.writerow(list(cluster_genes[k]))

    with open(os.path.join(high_level_dir, "deltas.csv"), "w", newline="") as file:
        writer = csv.writer(file)
        for d in deltas:
            writer.writerow([cluster_labels[d]] + list(deltas[d]))

    # with open(out_path + '/high_level_compress/counts.csv', 'w', newline='') as file:
    #    writer = csv.writer(file)
    #    for c in counts:
    #        writer.writerow(c)

    with open(os.path.join(high_level_dir, "counts.csv"), "w", newline="") as file:
        writer = csv.writer(file)
        for c in compressed_counts:
            writer.writerow(c)

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
        low_level_dir = os.path.join(out_path, "low_level_compress_count_extension")
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
    # low_level_compress('counts')
    low_level_compress("counts")


def high_level_decompress(in_path, matrix_path):
    gene_map = defaultdict(lambda: [None, set()])
    cluster_map = defaultdict(lambda: set())
    counts = []
    deltas_file = os.path.join(in_path, "high_level_compress_count_extension", "deltas.csv")
    with open(deltas_file, newline="") as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            cluster = int(row[0])
            genes = set(map(int, row[1:])) if len(row) > 1 else set()
            gene_map[i] = [cluster, genes]
            cluster_map[cluster].add(i)
            num_cells = i

    cluster_genes_file = os.path.join(
        in_path, "high_level_compress_count_extension", "cluster_genes.csv"
    )
    with open(cluster_genes_file, newline="") as f:
        reader = csv.reader(f)
        for cluster, row in enumerate(reader):
            genes = set(map(int, row))
            for cell in cluster_map[cluster]:
                gene_map[cell][1] = gene_map[cell][1].union(genes)

    counts_file = os.path.join(in_path, "high_level_compress_count_extension", "counts.csv")
    with open(counts_file, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            row_expanded = []
            for token in row:
                token = token.strip()
                if token == "":
                    continue
                try:
                    n = int(token)
                except ValueError:
                    # skip malformed tokens
                    continue
                if n < 0:
                    # negative tokens represent runs of 1s
                    row_expanded.extend([1] * (-n))
                else:
                    row_expanded.append(n)
            counts.append(row_expanded)

    # reconstruct into sparse COO to avoid dense memory
    original_sparse = mmread(matrix_path).tocsr()
    G, N = original_sparse.shape
    rows_list = []
    cols_list = []
    data_list = []
    # Diagnostic checks: ensure counts rows match expected non-zero gene counts per cell
    if len(counts) != (num_cells + 1):
        print(f"DEBUG: counts rows = {len(counts)}, expected = {num_cells + 1}")
    for cell_idx in range(0, num_cells + 1):
        gene_indices = sorted(gene_map[cell_idx][1])
        expected = len(gene_indices)
        got = len(counts[cell_idx]) if cell_idx < len(counts) else 0
        if expected != got:
            print(
                f"DEBUG MISMATCH cell {cell_idx}: expected {expected} values, got {got}"
            )
            print("  sample gene indices:", gene_indices[:20])
            print(
                "  sample counts:",
                counts[cell_idx][:40] if cell_idx < len(counts) else "MISSING",
            )
        for i, g in enumerate(gene_indices):
            val = counts[cell_idx][i] if i < got else 0
            if val != 0:
                rows_list.append(g)
                cols_list.append(cell_idx)
                data_list.append(val)

    reconstructed = coo_matrix(
        (data_list, (rows_list, cols_list)), shape=(G, N)
    ).tocsr()

    diff = (reconstructed != original_sparse).nnz
    passed = diff == 0
    print(
        "High Level Accuracy Check Passed"
        if passed
        else f"High Level Accuracy Check Failed (diff nnz={diff})"
    )
    return reconstructed, original_sparse


def low_level_decompress(in_path):
    def helper_decompress(target):
        tree_file = os.path.join(
            in_path, "low_level_compress_count_extension", f"{target}_huffman_tree"
        )
        with open(tree_file, "rb") as file:
            # Load the pickled object from the file
            root, length = pickle.load(file)

        decoded_lines = []
        count = 0
        enc_file = os.path.join(
            in_path, "low_level_compress_count_extension", f"huffman_encoded_{target}"
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
        hl_file = os.path.join(in_path, "high_level_compress_count_extension", f"{target}.csv")
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
