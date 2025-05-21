import os
import csv
import numpy as np

# Root directory to search for deltas.csv files
ROOT = os.path.dirname(os.path.abspath(__file__))

# List to store results for all deltas.csv files
results = []


def analyze_deltas_csv(filepath):
    """Analyze a deltas.csv file and return cluster statistics."""
    cluster_counts = {}
    with open(filepath, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            cluster = row[0]
            cluster_counts[cluster] = cluster_counts.get(cluster, 0) + 1
    sizes = list(cluster_counts.values())
    if not sizes:
        return None
    arr = np.array(sizes)
    return {
        "file": filepath,
        "num_clusters": len(cluster_counts),
        "avg": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=0)),
        "min": int(np.min(arr)),
        "max": int(np.max(arr)),
        "sizes": sizes,
    }


def find_deltas_csv_files(root):
    """Recursively find all deltas.csv files under root."""
    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            if fname == "deltas.csv":
                yield os.path.join(dirpath, fname)


def main():
    print("file,num_clusters,avg,std,min,max")
    for filepath in find_deltas_csv_files(ROOT):
        stats = analyze_deltas_csv(filepath)
        if stats:
            print(
                f"{os.path.relpath(stats['file'], ROOT)},{stats['num_clusters']},{stats['avg']:.2f},{stats['std']:.2f},{stats['min']},{stats['max']}"
            )


if __name__ == "__main__":
    main()
