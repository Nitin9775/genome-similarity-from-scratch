import time

from src.sketching import sketch_fasta_with_minimizers
from src.weighted_minimizers import weighted_minimizers


if __name__ == "__main__":
    genome = "data/ecoli.fasta"
    k = 15
    w = 10

    with open(genome) as f:
        sequence = "".join(
            line.strip()
            for line in f
            if not line.startswith(">")
        )

    print("Ordinary vs Weighted Minimizer Benchmark: E. coli")
    print(f"k={k}, w={w}")
    print()

    start = time.perf_counter()
    ordinary = list(
        sketch_fasta_with_minimizers(
            genome,
            k=k,
            w=w,
        ).hashes
    )
    ordinary_time = time.perf_counter() - start

    start = time.perf_counter()
    weighted = list(
        weighted_minimizers(
            sequence,
            k=k,
            w=w,
        )
    )
    weighted_time = time.perf_counter() - start

    print(f"Ordinary minimizers:  {len(ordinary):,}")
    print(f"Ordinary runtime:     {ordinary_time:.4f} seconds")
    print()

    print(f"Weighted minimizers:  {len(weighted):,}")
    print(f"Weighted runtime:     {weighted_time:.4f} seconds")