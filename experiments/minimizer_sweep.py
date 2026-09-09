import time

from src.sketching import kmerize, sketch_fasta_with_minimizers


if __name__ == "__main__":
    ecoli = "data/ecoli.fasta"
    shigella = "data/shigella.fasta"

    with open(ecoli) as f:
        sequence = "".join(
            line.strip()
            for line in f
            if not line.startswith(">")
        )

    total_kmers = sum(1 for _ in kmerize(sequence, 15))

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
        reduction = 100 * (1 - len(ecoli_sketch.hashes) / total_kmers)
        

        elapsed = time.perf_counter() - start

        print(f"w={w}")
        print(f"E. coli minimizers: {len(ecoli_sketch.hashes):,}")
        print(f"Shigella minimizers: {len(shigella_sketch.hashes):,}")
        print(f"Jaccard similarity: {jaccard}")
        print(f"E. coli reduction: {reduction:.2f}%")
        print(f"Jaccard change from w=5: {jaccard - 0.533537931030026:.6f}")
        print(f"Relative Jaccard change: {(jaccard - 0.533537931030026) / 0.533537931030026 * 100:.2f}%")
        print(f"Runtime: {elapsed:.4f} seconds")
        print()