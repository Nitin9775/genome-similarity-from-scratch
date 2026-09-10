import argparse
from email import parser
from src.sketching import sketch_fasta, sketch_fasta_with_minimizers
from src.ani_estimator import mash_distance


def main():
    parser = argparse.ArgumentParser(
        description="Estimate genome similarity using MinHash sketches."
    )

    parser.add_argument("genome1", help="Path to first FASTA file")
    parser.add_argument("genome2", help="Path to second FASTA file")
    parser.add_argument("--k", type=int, default=15, help="k-mer size")
    parser.add_argument("--size", type=int, default=1000, help="Sketch size")
    parser.add_argument(
        "--minimizers",
        action="store_true",
        help="Use minimizer sketches instead of MinHash",
    )
    parser.add_argument(
    "--w",
    type=int,
    default=10,
    help="Minimizer window size",
    )

    args = parser.parse_args()

    if args.minimizers:
        sketch1 = sketch_fasta_with_minimizers(
            args.genome1,
            k=args.k,
            w=args.w,
        )
        sketch2 = sketch_fasta_with_minimizers(
            args.genome2,
            k=args.k,
            w=args.w,
        )

    else:
        sketch1 = sketch_fasta(
            args.genome1,
            k=args.k,
            size=args.size,
        )
        sketch2 = sketch_fasta(
            args.genome2,
            k=args.k,
            size=args.size,
        )

    jaccard = sketch1.jaccard(sketch2)
    distance = mash_distance(jaccard, args.k)

    print(f"Jaccard similarity: {jaccard}")
    print(f"Mash distance: {distance}")


if __name__ == "__main__":
    main()