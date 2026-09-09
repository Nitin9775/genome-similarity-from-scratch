from src.sketching import exact_jaccard, sketch_fasta


if __name__ == "__main__":
    genome1 = "data/test.fasta"
    genome2 = "data/test_genome_2.fasta"

    k = 15

    exact = exact_jaccard(
        open(genome1).read().replace("\n", "").replace("\r", ""),
        open(genome2).read().replace("\n", "").replace("\r", ""),
        k,
    )

    print(f"Exact Jaccard: {exact}")

    for size in [10, 25, 50, 100, 500, 1000]:
        sketch1 = sketch_fasta(genome1, k=k, size=size)
        sketch2 = sketch_fasta(genome2, k=k, size=size)

        estimated = sketch1.jaccard(sketch2)
        error = abs(exact - estimated)

        print(
            f"Sketch size: {size:4d} | "
            f"MinHash: {estimated:.6f} | "
            f"Absolute error: {error:.6f}"
        )