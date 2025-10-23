import pickle
from collections import defaultdict
import csv
import numpy as np
from scipy.io import mmread
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

    # Extended compression of deltas
    if not os.path.isdir(out_path + '/high_level_compress2'):
        os.mkdir(out_path + '/high_level_compress2')
    copy_file(out_path + '/high_level_compress/cluster_genes.csv', out_path + '/high_level_compress2')
    copy_file(out_path + '/high_level_compress/counts.csv', out_path + '/high_level_compress2')

    # compress deltas further
    def compress_high_level_deltas(in_path):
        with open(in_path + '/high_level_compress/deltas.csv', newline='') as f:
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
        
        with open(in_path + '/high_level_compress2/deltas.csv', 'w', newline='') as file:
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

    compress_high_level_deltas(out_path)

    def low_level_compress(folder_num):
        for file_name in os.listdir(out_path + '/high_level_compress'+folder_num):
            if file_name != 'cluster_genes.csv':
                with open(out_path + '/high_level_compress'+folder_num+'/' + file_name, newline='') as f:
                    lines = f.readlines()

                    content = []
                    for l in lines:
                        content += list(map(str, l.replace('\r', '').replace('\n', '').split(',')))
                        content.append('\n')
                    content = content[:-1]

                    root = build_huffman_tree(collections.Counter(content))
                    huffman_codes = generate_huffman_codes(root)

                    # Pickle the array and save to a file
                    if not os.path.isdir(out_path + '/low_level_compress'+folder_num):
                        os.mkdir(out_path + '/low_level_compress'+folder_num)
                    filename = out_path + '/low_level_compress'+folder_num+'/' + file_name.replace('.csv', '') + '_huffman_tree'
                    with open(filename, 'wb') as file:
                        pickle.dump((root, len(content)), file)
                    file.close()
                    
                    with open(out_path + '/low_level_compress'+folder_num+'/huffman_encoded_'+file_name.replace('.csv', ''), 'wb') as file:
                        leftover = ''
                        for index, l in enumerate(lines):
                            elems = l.replace('\r', '').replace('\n', '').split(',')
                            elems = list(map(str, elems))
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
    
    low_level_compress('') # low-level compress original implementation
    low_level_compress('2') # low-level compress extended implementation

def high_level_decompress(in_path, matrix_path):
    gene_map = defaultdict(lambda: [None, set()])
    cluster_map = defaultdict(lambda: set())
    counts = []
    with open(in_path + '/deltas.csv', newline='') as f:
        reader = csv.reader(f)
        if in_path[-1] == '2':
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
        else:
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

def low_level_decompress(in_path, folder_num):
    def helper_decompress(target):
        with open(in_path + '/low_level_compress'+folder_num+'/' + target + '_huffman_tree', 'rb') as file:
            # Load the pickled object from the file
            root, length = pickle.load(file)
        file.close()

        decoded_lines = []
        count = 0
        with open(in_path + '/low_level_compress'+folder_num+'/huffman_encoded_' + target, 'rb') as file:
            decoded_line = []
            cur = root
            n = 0
            while byte := file.read(1):
                for i in range(7,-1,-1):
                    if count == length:
                        break
                    sym = cur.symbol
                    if sym is not None:
                        count += 1
                        if sym != '\n':
                            decoded_line.append(sym)
                        else:
                            decoded_lines.append(decoded_line)
                            decoded_line = []
                        cur = root
                    if (byte[0] >> i) & 1 == 0:
                        cur = cur.left
                    else:
                        cur = cur.right
            if decoded_line:
                decoded_lines.append(decoded_line)
        file.close()

        result = 'ENCODED CORRECTLY!'
        lines = ''
        with open(in_path + '/high_level_compress'+folder_num+'/' + target + '.csv', newline='') as f:
            lines = f.readlines()
        for i,l in enumerate(lines):
            l = list(l.replace('\r', '').replace('\n', '').split(','))
            if (l != decoded_lines[i]):
                print(l)
                print(decoded_lines[i])
                result = 'ENCODED WRONG!!'
                break
        print(target + ' huffman : ' + result)

    helper_decompress('deltas')
    helper_decompress('counts')