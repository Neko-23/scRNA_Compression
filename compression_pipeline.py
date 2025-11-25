import os
from neural import assign_clusters
from scipy.io import mmread
from scipy.sparse import save_npz
import numpy as np
from compressor import compress
from compressor_delta_extension import compress as extra_delta_compress
from compressor_counts_extension import compress as extra_count_compress
from compressor_counts_and_delta_extension import compress as extra_count_and_delta_compress
from visualization import create_bar_charts, create_pie_chart
import sys

def main(sample, k):
    print(sample)
    print(k)
    label_path = 'compressed/neural/sample'+str(sample)+'/k'+str(k)+'/cluster_labels.pkl'
    matrix_path = 'data/sample'+str(sample)+'/matrix.mtx'
    out_path = 'compressed/neural/sample'+str(sample)+'/k'+str(k)

    os.makedirs('compressed/neural/sample'+str(sample), exist_ok=True)
    os.makedirs(out_path, exist_ok=True)

    assign_clusters(sample, k)

    sparse_matrix = mmread(matrix_path)
    csr = sparse_matrix.tocsr()
    save_npz('data/sample'+str(sample)+'/csr.npz', csr)
    output_prefix = 'data/sample'+str(sample)+'/csr_nozip/csr'
    if not os.path.isdir('data/sample'+str(sample)+'/csr_nozip'):
        os.mkdir('data/sample'+str(sample)+'/csr_nozip')
    np.save(f"{output_prefix}_data.npy", csr.data)
    np.save(f"{output_prefix}_indices.npy", csr.indices)
    np.save(f"{output_prefix}_indptr.npy", csr.indptr)
    np.save(f"{output_prefix}_shape.npy", np.array(csr.shape))

    csc = sparse_matrix.tocsc()
    save_npz('data/sample'+str(sample)+'/csc.npz', csc)
    output_prefix = 'data/sample'+str(sample)+'/csc_nozip/csc'
    if not os.path.isdir('data/sample'+str(sample)+'/csc_nozip'):
        os.mkdir('data/sample'+str(sample)+'/csc_nozip')
    np.save(f"{output_prefix}_data.npy", csc.data)
    np.save(f"{output_prefix}_indices.npy", csc.indices)
    np.save(f"{output_prefix}_indptr.npy", csc.indptr)
    np.save(f"{output_prefix}_shape.npy", np.array(csc.shape))

    compress(label_path, matrix_path, out_path)
    extra_delta_compress(label_path, matrix_path, out_path)
    extra_count_compress(label_path, matrix_path, out_path)
    extra_count_and_delta_compress(label_path, matrix_path, out_path)

    import shutil

    # copy cluster_genes.csv from high level to low level
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
        
    copy_file(out_path + '/high_level_compress/cluster_genes.csv', out_path + '/low_level_compress')
    copy_file(out_path + '/high_level_compress_count_extension/cluster_genes.csv', out_path + '/low_level_compress_count_extension')
    copy_file(out_path + '/high_level_compress_delta_extension/cluster_genes.csv', out_path + '/low_level_compress_delta_extension')
    copy_file(out_path + '/high_level_compress_count_and_delta_extension/cluster_genes.csv', out_path + '/low_level_compress_count_and_delta_extension')

    # gzip high_level_compress dir and low_level_compress dir
    import subprocess
    subprocess.run(["tar", "-czf", out_path + '/high_level_compress.tar.gz', out_path + '/high_level_compress'])
    subprocess.run(["tar", "-czf", out_path + '/low_level_compress.tar.gz', out_path + '/low_level_compress'])
    subprocess.run(["tar", "-czf", out_path + '/high_level_compress_delta_extension.tar.gz', out_path + '/high_level_compress_delta_extension'])
    subprocess.run(["tar", "-czf", out_path + '/low_level_compress_delta_extension.tar.gz', out_path + '/low_level_compress_delta_extension'])
    subprocess.run(["tar", "-czf", out_path + '/high_level_compress_count_extension.tar.gz', out_path + '/high_level_compress_count_extension'])
    subprocess.run(["tar", "-czf", out_path + '/low_level_compress_count_extension.tar.gz', out_path + '/low_level_compress_count_extension'])
    subprocess.run(["tar", "-czf", out_path + '/high_level_compress_count_and_delta_extension.tar.gz', out_path + '/high_level_compress_count_and_delta_extension'])
    subprocess.run(["tar", "-czf", out_path + '/low_level_compress_count_and_delta_extension.tar.gz', out_path + '/low_level_compress_count_and_delta_extension'])
    subprocess.run(["tar", "-czf", 'data/sample'+str(sample)+'/matrix.tar.gz', 'data/sample'+str(sample)+'/matrix.mtx'])

    from compressor import high_level_decompress, low_level_decompress
    from compressor_delta_extension import high_level_decompress as hld_delta
    from compressor_delta_extension import low_level_decompress as lld_delta
    from compressor_counts_extension import high_level_decompress as hld_counts
    from compressor_counts_extension import low_level_decompress as lld_counts
    from compressor_counts_and_delta_extension import high_level_decompress as hld_counts_and_delta
    from compressor_counts_and_delta_extension import low_level_decompress as lld_counts_and_delta

    in_path = 'compressed/neural/sample'+str(sample)+'/k'+str(k)
    high_level_decompress(in_path + '/high_level_compress', matrix_path)
    low_level_decompress(in_path)
    hld_delta(in_path, matrix_path)
    lld_delta(in_path)
    hld_counts(in_path, matrix_path)
    lld_counts(in_path)
    hld_counts_and_delta(in_path, matrix_path)
    lld_counts_and_delta(in_path)

    create_bar_charts(out_path, sample, k)
    create_pie_chart(out_path, sample, k)

if __name__ == "__main__":
    sample = int(sys.argv[1])
    k = int(sys.argv[2])

    main(sample, k)