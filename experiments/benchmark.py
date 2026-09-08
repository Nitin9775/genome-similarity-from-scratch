import time

from src.sketching import sketch_fasta
from src.ani_estimator import mash_distance

def benchmark_genome_pair(genome1, genome2, k=15, size=1000):
    start = time.perf_counter()

    sketch1 = sketch_fasta(genome1, k=k, size=size)
    sketch2 = sketch_fasta(genome2, k=k, size=size)

    jaccard = sketch1.jaccard(sketch2)
    distance = mash_distance(jaccard, k)

    elapsed = time.perf_counter() - start

    return jaccard, distance, elapsed

if __name__ == "__main__":
    genome1 = "data/ecoli.fasta"
    genome2 = "data/shigella.fasta"

    jaccard, distance, elapsed = benchmark_genome_pair(
        genome1,
        genome2,
        k=15,
        size=1000,
    )

    print(f"Jaccard similarity: {jaccard}")
    print(f"Mash distance: {distance}")
    print(f"Runtime: {elapsed:.4f} seconds")