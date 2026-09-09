import time

from src.sketching import sketch_fasta_with_minimizers


if __name__ == "__main__":
    ecoli = "data/ecoli.fasta"
    shigella = "data/shigella.fasta"

    for w in [5, 10, 20, 50, 100]:
        start = time.perf_counter()

        ecoli_sketch = sketch_fasta_with_minimizers(
            ecoli,
            k=15,
            w=w,
        )

        shigella_sketch = sketch_fasta_with_minimizers(
            shigella,
            k=15,
            w=w,
        )

        jaccard = ecoli_sketch.jaccard(shigella_sketch)

        elapsed = time.perf_counter() - start

        print(f"w={w}")
        print(f"E. coli minimizers: {len(ecoli_sketch.hashes):,}")
        print(f"Shigella minimizers: {len(shigella_sketch.hashes):,}")
        print(f"Jaccard similarity: {jaccard}")
        print(f"Runtime: {elapsed:.4f} seconds")
        print()