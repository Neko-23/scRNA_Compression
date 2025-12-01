import pickle
from collections import defaultdict
import csv
import numpy as np
from scipy.io import mmread
from scipy.sparse import coo_matrix, csr_matrix
import heapq
import collections
import os
import time
import psutil
import math

# fpgrowth libraries
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


def compress(cluster_assignments_file, matrix_path, out_path):
    process = psutil.Process(os.getpid())
    start = int(process.memory_info().rss)
    start_time = time.perf_counter()

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

    mid_time = time.perf_counter()
    mid_rss = int(process.memory_info().rss)
    runlentoks = np.concatenate(compressed_counts)
    print(f'total runlentok: {np.sum(runlentoks)}')
    print(f'median runlentok: {np.median(runlentoks)}')
    print(f'75th percentile runlentok: {np.percentile(runlentoks, 75)}')
    print(f'95th percentile runlentok: {np.percentile(runlentoks, 95)}')
    mid_time2 = time.perf_counter()
    mid_rss2 = int(process.memory_info().rss)

    high_level_dir = os.path.join(out_path, "high_level_compress_count_and_delta_extension")
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

    with open(os.path.join(high_level_dir, "counts.csv"), "w", newline="") as file:
        writer = csv.writer(file)
        for c in compressed_counts:
            writer.writerow(c)
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time - (mid_time2 - mid_time)
    current = int(process.memory_info().rss)
    elapsed_rss = int(current)-int(start)
    elapsed_rss -= int(mid_rss2-mid_rss)
    print(f'high-level rss (MB): {elapsed_rss/1e6}')
    print(f"high-level elapsed wall time: {elapsed_time:.6f} seconds")
    
    # compress deltas further
    def compress_high_level_deltas():
        with open(os.path.join(high_level_dir, 'deltas.csv'), newline='') as f:
            reader = csv.reader(f)
            dataset = []
            cluster_labels = []
            for i, row in enumerate(reader):
                dataset.append(row[1:])
                cluster_labels.append(row[0])
            f.close()

        chosens = []
        flags = []
        support = 0.95
        iter = 1
        while len(chosens) < 2000:
            te = TransactionEncoder()
            te_ary = te.fit(dataset).transform(dataset)
            df = pd.DataFrame(te_ary, columns=te.columns_)
            
            if support < 1:
                df_iter = df.sample(n = 2, random_state=iter)
            else:
                df_iter = df

            s = fpgrowth(df_iter, min_support=support, use_colnames=True, max_len = 2)

            mlength = 1
            chosen = []
            choice = None
            for ele in s['itemsets']:
                if len(ele) > mlength:
                    mlength = len(ele)
                    choice = set(ele)
            if choice == None:
                support -= 0.01
                support = round(support, 2)
                iter += 1
                continue
            chosen.append(choice)
            temp_set = choice
            for ele in s['itemsets']:
                ele_set = set(ele)
                if len(ele) > 1 and not temp_set.intersection(ele_set):
                    temp_set = temp_set.union(ele_set)
                    chosen.append(ele_set)
            if len(chosen) < 50:
                support -= 0.01
                support = round(support, 2)
            chosens = chosens + chosen

            deltas = [set(x) for x in dataset]
            flag = []
            for c in chosen:
                f = set()
                for i, d in enumerate(deltas):
                    if c.intersection(d) == c:
                        deltas[i] = (d - c)
                        f.add(i)
                flag.append(f)
            flags += flag
            print('support: ', support)
            print('chosens: ', len(chosens))
            iter += 1
            dataset = [list(x) for x in deltas]
        
        with open(os.path.join(high_level_dir, 'deltas.csv'), 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow([len(chosens)]) # number of common item sets
            for chosen in chosens:
                writer.writerow(list(chosen)) # write current commom item set
            for i, d in enumerate(deltas):
                row = [cluster_labels[i]]
                for j, flag in enumerate(flags):
                    if i in flag:
                        row += ['#'+str(j)]
                row += list(d)
                writer.writerow(row)
            file.close()

    process = psutil.Process(os.getpid())
    start = int(process.memory_info().rss)
    start_time = time.perf_counter()
    compress_high_level_deltas()
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    current = int(process.memory_info().rss)
    print(f'fpgrowth rss (MB): {(current-start)/1e6}')
    print(f"fpgrowth elapsed wall time: {elapsed_time:.6f} seconds")

    def low_level_compress(target):
        src_file = os.path.join(high_level_dir, f"{target}.csv")
        with open(src_file, newline="") as f:
            lines = f.readlines()

        content = []
        for l in lines:
            # tolerate stray whitespace and empty tokens
            content += list(map(str, l.replace("\r", "").replace("\n", "").split(",")))
            content.append("\n")
        if content:
            content = content[:-1]

        freq_dict = collections.Counter(content)
        root = build_huffman_tree(freq_dict)
        huffman_codes = generate_huffman_codes(root)
        
        mid_time = time.perf_counter()
        mid_rss = int(process.memory_info().rss)
        print(f'{target} # distinct symbols: {len(huffman_codes.keys())}')
        print(f'{target} avg code len: {np.mean([len(x) for x in huffman_codes.values()])}')
        # calculate entropy
        total = sum(freq_dict.values())
        entropy = 0
        for v in freq_dict.values():
            prob = v/total
            entropy += (-prob*math.log2(prob))
        print(f'{target} entropy: {entropy}')
        mid_time2 = time.perf_counter()
        mid_rss2 = int(process.memory_info().rss)
            
        # Pickle the array and save to a file
        low_level_dir = os.path.join(out_path, "low_level_compress_count_and_delta_extension")
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
                    elems = list(map(str, elems))
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
        
        return (mid_time2-mid_time, mid_rss2-mid_rss)

    process = psutil.Process(os.getpid())
    start_time = time.perf_counter()
    start = int(process.memory_info().rss)
    mid_time_duration, mid_rss_chunk = low_level_compress("deltas")
    mid_time_duration2, mid_rss_chunk2 = low_level_compress("counts")
    current = int(process.memory_info().rss)
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time - mid_time_duration - mid_time_duration2
    elapsed_rss = current-start
    print(f'huffman rss (MB): {elapsed_rss/1e6}')
    print(f"huffman elapsed wall time: {elapsed_time:.6f} seconds")


def high_level_decompress(in_path, matrix_path):
    gene_map = defaultdict(lambda: [None, set()])
    cluster_map = defaultdict(lambda: set())
    counts = []
    deltas_file = os.path.join(in_path, "high_level_compress_count_and_delta_extension", "deltas.csv")
    with open(deltas_file, newline='') as f:
        reader = csv.reader(f)
        num_item_sets = 0
        item_sets = []
        deltas = []
        for i, row in enumerate(reader):
            if i == 0:
                num_item_sets = int(row[0])
            elif i <= num_item_sets:
                item_sets.append((int(row[0]), int(row[1])))
            else:
                delta = []
                for c in row:
                    if c[0] != '#':
                        delta.append(int(c))
                    else:
                        a, b = item_sets[int(c[1:])]
                        delta.append(a)
                        delta.append(b)
                deltas.append(delta)
        
        for i, d in enumerate(deltas):
            gene_map[i] = [d[0], set(d[1:])]
            cluster_map[d[0]].add(i)
            num_cells = i            
        f.close()

    cluster_genes_file = os.path.join(
        in_path, "high_level_compress_count_and_delta_extension", "cluster_genes.csv"
    )
    with open(cluster_genes_file, newline="") as f:
        reader = csv.reader(f)
        for cluster, row in enumerate(reader):
            genes = set(map(int, row))
            for cell in cluster_map[cluster]:
                gene_map[cell][1] = gene_map[cell][1].union(genes)

    counts_file = os.path.join(in_path, "high_level_compress_count_and_delta_extension", "counts.csv")
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
            in_path, "low_level_compress_count_and_delta_extension", f"{target}_huffman_tree"
        )
        with open(tree_file, "rb") as file:
            # Load the pickled object from the file
            root, length = pickle.load(file)

        decoded_lines = []
        count = 0
        enc_file = os.path.join(
            in_path, "low_level_compress_count_and_delta_extension", f"huffman_encoded_{target}"
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
        lines = ''
        hl_file = os.path.join(in_path, "high_level_compress_count_and_delta_extension", f"{target}.csv")
        with open(hl_file, newline="") as f:
            lines = f.readlines()
        for i, l in enumerate(lines):
            l = list(l.replace('\r', '').replace('\n', '').split(','))
            if (l != decoded_lines[i]):
                print(l)
                print(decoded_lines[i])
                result = 'ENCODED WRONG!!'
                break
        print(target + " huffman : " + result)

    helper_decompress("deltas")
    helper_decompress("counts")
