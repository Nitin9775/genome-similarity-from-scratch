import argparse

from src.sketching import sketch_fasta
from src.ani_estimator import mash_distance


def main():
    parser = argparse.ArgumentParser(
        description="Estimate genome similarity using MinHash sketches."
    )

    parser.add_argument("genome1", help="Path to first FASTA file")
    parser.add_argument("genome2", help="Path to second FASTA file")
    parser.add_argument("--k", type=int, default=15, help="k-mer size")
    parser.add_argument("--size", type=int, default=1000, help="Sketch size")

    args = parser.parse_args()

    sketch1 = sketch_fasta(args.genome1, k=args.k, size=args.size)
    sketch2 = sketch_fasta(args.genome2, k=args.k, size=args.size)

    jaccard = sketch1.jaccard(sketch2)
    distance = mash_distance(jaccard, args.k)

    print(f"Jaccard similarity: {jaccard}")
    print(f"Mash distance: {distance}")


if __name__ == "__main__":
    main()