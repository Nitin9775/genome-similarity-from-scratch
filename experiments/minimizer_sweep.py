import time

from src.sketching import sketch_fasta_with_minimizers


if __name__ == "__main__":
    genome = "data/ecoli.fasta"

    for w in [5, 10, 20, 50, 100]:
        start = time.perf_counter()

        sketch = sketch_fasta_with_minimizers(
            genome,
            k=15,
            w=w,
        )

        elapsed = time.perf_counter() - start

        print(f"w={w}")
        print(f"Minimizers: {len(sketch.hashes)}")
        print(f"Runtime: {elapsed:.4f} seconds")
        print()