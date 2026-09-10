from src.sketching import exact_jaccard, sketch_fasta
from src.sketching import sketch_fasta_with_minimizers


if __name__ == "__main__":
    genome1 = "data/ecoli.fasta"
    genome2 = "data/shigella.fasta"

    k = 15
    sketch_size = 1000
    w = 10

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

    minhash1 = sketch_fasta(genome1, k=k, size=sketch_size)
    minhash2 = sketch_fasta(genome2, k=k, size=sketch_size)
    minhash_jaccard = minhash1.jaccard(minhash2)

    minimizer1 = sketch_fasta_with_minimizers(
        genome1,
        k=k,
        w=w,
    )
    minimizer2 = sketch_fasta_with_minimizers(
        genome2,
        k=k,
        w=w,
    )
    minimizer_jaccard = minimizer1.jaccard(minimizer2)

    print("MinHash vs Minimizer Accuracy")
    print(f"k={k}, sketch_size={sketch_size}, w={w}")
    print()

    print(f"Exact Jaccard:       {exact:.6f}")
    print(f"MinHash Jaccard:     {minhash_jaccard:.6f}")
    print(f"MinHash error:       {abs(exact - minhash_jaccard):.6f}")
    print()
    print(f"Minimizer Jaccard:   {minimizer_jaccard:.6f}")
    print(f"Minimizer error:     {abs(exact - minimizer_jaccard):.6f}")