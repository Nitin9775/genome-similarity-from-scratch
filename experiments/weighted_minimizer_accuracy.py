from src.sketching import exact_jaccard, sketch_fasta_with_minimizers
from src.weighted_minimizers import weighted_minimizers


def sequence_from_fasta(filepath):
    with open(filepath) as f:
        return "".join(
            line.strip()
            for line in f
            if not line.startswith(">")
        )


if __name__ == "__main__":
    genome1 = "data/ecoli.fasta"
    genome2 = "data/shigella.fasta"

    k = 15
    w = 10

    sequence1 = sequence_from_fasta(genome1)
    sequence2 = sequence_from_fasta(genome2)

    exact = exact_jaccard(sequence1, sequence2, k)

    ordinary1 = sketch_fasta_with_minimizers(
        genome1,
        k=k,
        w=w,
    )
    ordinary2 = sketch_fasta_with_minimizers(
        genome2,
        k=k,
        w=w,
    )

    weighted1 = set(weighted_minimizers(sequence1, k=k, w=w))
    weighted2 = set(weighted_minimizers(sequence2, k=k, w=w))

    ordinary_jaccard = ordinary1.jaccard(ordinary2)

    weighted_jaccard = (
        len(weighted1 & weighted2) /
        len(weighted1 | weighted2)
    )

    print("Ordinary vs Weighted Minimizer Accuracy")
    print(f"k={k}, w={w}")
    print()
    print(f"Exact Jaccard:              {exact:.6f}")
    print()
    print(f"Ordinary minimizer Jaccard: {ordinary_jaccard:.6f}")
    print(f"Ordinary error:             {abs(exact - ordinary_jaccard):.6f}")
    print()
    print(f"Weighted minimizer Jaccard: {weighted_jaccard:.6f}")
    print(f"Weighted error:             {abs(exact - weighted_jaccard):.6f}")