from src.sketching import exact_jaccard, sketch_fasta_with_minimizers


if __name__ == "__main__":
    genome1 = "data/ecoli.fasta"
    genome2 = "data/shigella.fasta"

    k = 15

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

    print("Minimizer Accuracy Sweep")
    print(f"k: {k}")
    print(f"Exact Jaccard: {exact:.6f}")
    print()

    for w in [5, 10, 20, 50, 100]:
        sketch1 = sketch_fasta_with_minimizers(
            genome1,
            k=k,
            w=w,
        )
        sketch2 = sketch_fasta_with_minimizers(
            genome2,
            k=k,
            w=w,
        )

        estimated = sketch1.jaccard(sketch2)
        error = abs(exact - estimated)

        print(f"Window size: {w:3d}")
        print(f"  E. coli minimizers:   {len(sketch1.hashes):,}")
        print(f"  Shigella minimizers:  {len(sketch2.hashes):,}")
        print(f"  MinHash Jaccard:       {estimated:.6f}")
        print(f"  Absolute error:        {error:.6f}")
        print()