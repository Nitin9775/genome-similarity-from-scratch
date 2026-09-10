import time

from src.sketching import sketch_fasta


if __name__ == "__main__":
    genome = "data/ecoli.fasta"

    print("MinHash Runtime vs k-mer Size: E. coli")
    print()

    for k in [15, 21, 31, 41]:
        start = time.perf_counter()

        sketch = sketch_fasta(
            genome,
            k=k,
            size=1000,
        )

        elapsed = time.perf_counter() - start

        print(
            f"k={k:2d} | "
            f"Stored hashes: {len(sketch.hashes):4d} | "
            f"Runtime: {elapsed:.4f} seconds"
        )