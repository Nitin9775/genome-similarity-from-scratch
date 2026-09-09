## Current Implementation

![Tests](https://img.shields.io/badge/tests-66%20passed-brightgreen)
[![Tests](https://github.com/Nitin9775/genome-similarity-from-scratch/actions/workflows/tests.yml/badge.svg)](https://github.com/Nitin9775/genome-similarity-from-scratch/actions/workflows/tests.yml)

The project currently supports:

- Canonical k-mer generation
- FASTA file reading
- Deterministic hashing
- Bottom-k MinHash sketches
- Jaccard similarity estimation
- Mash distance calculation
- Streaming minimizer generation
- Minimizer-based sketching

### Example

```python
from src.sketching import sketch_fasta
from src.ani_estimator import mash_distance

sketch1 = sketch_fasta("data/test.fasta", k=15, size=1000)
sketch2 = sketch_fasta("data/test_genome_2.fasta", k=15, size=1000)

jaccard = sketch1.jaccard(sketch2)
distance = mash_distance(jaccard, k=15)

print("Jaccard similarity:", jaccard)
print("Mash distance:", distance)


### Canonical k-mer Generation

The implementation generates canonical k-mers by considering both a
k-mer and its reverse complement and selecting the lexicographically
smaller sequence.

Invalid DNA k-mers containing characters other than A, C, G, and T
are skipped.

### MinHash Sketching

A bottom-k MinHash sketch is implemented from scratch using deterministic
MD5 hashing.

The sketch stores the smallest hash values up to a configurable sketch
size. Jaccard similarity is estimated from the resulting sketches.

### Mash Distance

Mash distance is calculated from the estimated Jaccard similarity:

D = -(1/k) ln(2J/(1+J))

where J is the Jaccard similarity and k is the k-mer size.

### FASTA Processing

The project includes a FASTA reader that extracts and combines sequence
lines while ignoring header lines.

The command-line interface can compare two genome FASTA files and report
their estimated Jaccard similarity and Mash distance.

Example:

```text
python -m src.cli data/ecoli.fasta data/shigella.fasta --k 15 --size 1000