import time

from src.sketching import sketch_fasta


if __name__ == "__main__":
    genome = "data/ecoli.fasta"

    print("MinHash benchmark: E. coli")
    print()

    for size in [100, 500, 1000, 5000]:
        start = time.perf_counter()

        sketch = sketch_fasta(
            genome,
            k=15,
            size=size,
        )

        elapsed = time.perf_counter() - start

        print(
            f"Sketch size: {size:5d} | "
            f"Stored hashes: {len(sketch.hashes):5d} | "
            f"Runtime: {elapsed:.4f} seconds"
        )