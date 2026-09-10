from src.sketching import exact_jaccard, sketch_sequence


if __name__ == "__main__":
    sequence1 = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTA" * 100
    sequence2 = "ATGCGTACGTAGCTAGCTAGATAGCTAGCTAGCTA" * 100

    k = 15

    exact = exact_jaccard(sequence1, sequence2, k)

    print(f"Sequence length: {len(sequence1):,}")
    print(f"k: {k}")
    print(f"Exact Jaccard: {exact:.6f}")
    print()

    for size in [50, 100, 500, 1000]:
        sketch1 = sketch_sequence(sequence1, k=k, size=size)
        sketch2 = sketch_sequence(sequence2, k=k, size=size)

        estimated = sketch1.jaccard(sketch2)
        error = abs(exact - estimated)

        print(
            f"Sketch size: {size:4d} | "
            f"MinHash: {estimated:.6f} | "
            f"Absolute error: {error:.6f}"
        )