import time

from src.sketching import sketch_fasta
from src.sketching import sketch_fasta_with_minimizers


if __name__ == "__main__":
    genome = "data/ecoli.fasta"
    k = 15
    sketch_size = 1000
    w = 10

    print("MinHash vs Minimizer Benchmark: E. coli")
    print(f"k={k}, sketch_size={sketch_size}, w={w}")
    print()

    start = time.perf_counter()
    minhash = sketch_fasta(genome, k=k, size=sketch_size)
    minhash_time = time.perf_counter() - start

    start = time.perf_counter()
    minimizer_sketch = sketch_fasta_with_minimizers(
        genome,
        k=k,
        w=w,
    )
    minimizer_time = time.perf_counter() - start

    print(f"MinHash stored hashes:     {len(minhash.hashes)}")
    print(f"MinHash runtime:            {minhash_time:.4f} seconds")
    print()

    print(f"Minimizer stored hashes:    {len(minimizer_sketch.hashes)}")
    print(f"Minimizer runtime:          {minimizer_time:.4f} seconds")
    print()

    print(
        f"Runtime ratio (MinHash / Minimizer): "
        f"{minhash_time / minimizer_time:.2f}x"
    )