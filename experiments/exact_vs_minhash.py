from src.sketching import exact_jaccard, sketch_fasta


if __name__ == "__main__":
    genome1 = "data/test.fasta"
    genome2 = "data/test_genome_2.fasta"

    k = 15
    sketch_size = 1000

    with open(genome1) as f:
        sequence1 = "".join(
            line.strip()
            for line in f
            if not line.startswith(">")
        )

    with open(genome2) as f:
        sequence2 = "".join(
            line.strip()
            for line in f
            if not line.startswith(">")
        )

    exact = exact_jaccard(sequence1, sequence2, k)

    sketch1 = sketch_fasta(genome1, k=k, size=sketch_size)
    sketch2 = sketch_fasta(genome2, k=k, size=sketch_size)

    estimated = sketch1.jaccard(sketch2)

    print(f"Exact Jaccard: {exact}")
    print(f"MinHash Jaccard: {estimated}")
    print(f"Absolute error: {abs(exact - estimated)}")