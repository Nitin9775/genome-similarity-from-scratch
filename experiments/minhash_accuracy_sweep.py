from src.sketching import exact_jaccard, sketch_sequence


def mutate_sequence(sequence, step):
    sequence = list(sequence)

    for i in range(0, len(sequence), step):
        sequence[i] = {
            "A": "C",
            "C": "G",
            "G": "T",
            "T": "A",
        }[sequence[i]]

    return "".join(sequence)


if __name__ == "__main__":
    base = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTA" * 100

    k = 15

    print("MinHash Accuracy Sweep")
    print(f"Sequence length: {len(base):,}")
    print(f"k: {k}")
    print()

    for step in [10, 20, 50, 100]:
        sequence2 = mutate_sequence(base, step)

        exact = exact_jaccard(base, sequence2, k)

        print(f"Mutation step: {step}")
        print(f"Exact Jaccard: {exact:.6f}")

        for size in [50, 100, 500, 1000]:
            sketch1 = sketch_sequence(base, k=k, size=size)
            sketch2 = sketch_sequence(sequence2, k=k, size=size)

            estimated = sketch1.jaccard(sketch2)
            error = abs(exact - estimated)

            print(
                f"  Sketch size: {size:4d} | "
                f"MinHash: {estimated:.6f} | "
                f"Absolute error: {error:.6f}"
            )

        print()