# Genome Similarity From Scratch

A research-oriented implementation of genome similarity estimation using
k-mers, MinHash sketches, Mash distance, and streaming minimizers.

The goal of this project is to understand and implement the core
algorithms behind genome similarity estimation rather than simply
relying on existing tools.

---

## Current Implementation

![Tests](https://img.shields.io/badge/tests-68%20passed-brightgreen)

[![Tests](https://github.com/Nitin9775/genome-similarity-from-scratch/actions/workflows/tests.yml/badge.svg)](https://github.com/Nitin9775/genome-similarity-from-scratch/actions/workflows/tests.yml)

The project currently supports:

- Canonical k-mer generation
- Reverse-complement handling
- FASTA file reading
- Deterministic hashing
- Bottom-k MinHash sketches
- Jaccard similarity estimation
- Mash distance calculation
- Streaming minimizer generation
- Minimizer-based sketching
- Parameter-sweep experiments
- Automated testing

---

## Methods

### Canonical k-mers

A DNA sequence is divided into overlapping k-mers.

For each k-mer, both the forward sequence and its reverse complement
are considered. The lexicographically smaller sequence is selected as
the canonical representation.

This prevents a sequence and its reverse complement from being treated
as different k-mers.

Invalid k-mers containing characters other than `A`, `C`, `G`, and `T`
are skipped.

### MinHash Sketching

The project implements a bottom-k MinHash sketch from scratch.

Each canonical k-mer is converted into a deterministic hash value.

Only the smallest hash values up to a configurable sketch size are
retained.

The resulting sketches are used to estimate Jaccard similarity between
two genomes.

The implementation uses deterministic MD5 hashing so that results remain
reproducible across Python processes.

### Jaccard Similarity

For two sets of k-mers, Jaccard similarity is defined as:

```text
J(A,B) = |A ∩ B| / |A ∪ B|
```

The project estimates this similarity from MinHash or minimizer sketches.

A value closer to `1` indicates greater similarity, while a value closer
to `0` indicates lower similarity.

### Mash Distance

Mash distance is calculated from the estimated Jaccard similarity:

```text
D = -(1/k) ln(2J/(1+J))
```

where:

* `D` = Mash distance
* `J` = Jaccard similarity
* `k` = k-mer size

---

## FASTA Processing

The project includes a FASTA reader that:

* Ignores FASTA header lines
* Handles multiple sequence lines
* Combines sequence lines into a single sequence
* Supports the genome FASTA files used in the experiments

---

## Command-Line Interface

Two genome FASTA files can be compared directly from the command line.

### Example

```bash
python -m src.cli data/ecoli.fasta data/shigella.fasta --k 15 --size 1000
```

### Example Output

```text
Jaccard similarity: 0.539
Mash distance: 0.02373502549132912
```

---

## Python Usage

The core functionality can also be used directly from Python:

```python
from src.sketching import sketch_fasta
from src.ani_estimator import mash_distance

sketch1 = sketch_fasta(
    "data/test.fasta",
    k=15,
    size=1000
)

sketch2 = sketch_fasta(
    "data/test_genome_2.fasta",
    k=15,
    size=1000
)

jaccard = sketch1.jaccard(sketch2)
distance = mash_distance(jaccard, k=15)

print("Jaccard similarity:", jaccard)
print("Mash distance:", distance)
```

---

## Minimizer-based Sketching

The project also implements streaming minimizer selection from scratch.

For each window of consecutive k-mers, the minimum canonical k-mer hash
is selected as the minimizer.

The implementation uses a deque-based streaming algorithm rather than
materializing the complete k-mer hash list.

Deterministic MD5 hashing is used for reproducibility.

---

## Minimizer Window Experiment

The effect of minimizer window size was evaluated using the E. coli K-12
and Shigella genomes with:

```text
k = 15
```

### E. coli Minimizer Reduction

| Window Size | E. coli Minimizers | Reduction | Retention |
| ----------: | -----------------: | --------: | --------: |
|           5 |          1,494,480 |    67.80% |    32.20% |
|          10 |            816,575 |    82.41% |    17.59% |
|          20 |            428,206 |    90.77% |     9.23% |
|          50 |            176,237 |    96.20% |     3.80% |
|         100 |             89,107 |    98.08% |     1.92% |

Increasing the minimizer window size substantially reduces the number of
selected minimizers.

At `w=100`, 98.08% of the canonical E. coli k-mers are removed while
1.92% are retained.

---

## E. coli vs Shigella Similarity

The effect of minimizer window size on observed similarity was also
measured.

| Window Size | E. coli Minimizers | Shigella Minimizers | Jaccard Similarity |
| ----------: | -----------------: | ------------------: | -----------------: |
|           5 |          1,494,480 |           1,436,830 |             0.5335 |
|          10 |            816,575 |             784,905 |             0.5278 |
|          20 |            428,206 |             412,133 |             0.5173 |
|          50 |            176,237 |             169,840 |             0.5023 |
|         100 |             89,107 |              85,667 |             0.4943 |

The results demonstrate a trade-off between sequence reduction and the
similarity signal retained by the sketch.

As the window size increases:

* The number of minimizers decreases
* The reduction in sequence representation increases
* The observed Jaccard similarity decreases

For this experiment:

```text
w = 5
Jaccard = 0.5335

w = 100
Jaccard = 0.4943
```

The relative Jaccard change from `w=5` to `w=100` is approximately:

```text
-7.36%
```

---

## MinHash vs Minimizer Comparison

For E. coli and Shigella with `k=15`:

```text
MinHash Jaccard similarity:
0.539

Minimizer Jaccard similarity (w=10):
0.5277813969464745
```

The difference in observed Jaccard similarity is:

```text
0.539 - 0.5277813969464745
= 0.0112186030535255
```

This provides a direct comparison between similarity estimated using
bottom-k MinHash sketches and similarity estimated from minimizer-based
sketches.

---

## Streaming Implementation

The minimizer implementation processes the sequence using a streaming
deque-based algorithm.

For the E. coli genome with:

```text
k = 15
w = 10
```

the recorded runtime comparison was:

```text
Full MinHash sketch:
57.17 seconds

Streaming minimizer generation:
14.68 seconds
```

### Performance

```text
Speedup:
3.90x

Runtime reduction:
74.33%
```

These measurements demonstrate the computational benefit of streaming
minimizer generation for this implementation.

---

## Parameter Experiments

The repository contains experiments investigating:

* k-mer size
* Minimizer window size
* Number of selected minimizers
* Minimizer reduction
* Minimizer retention
* Jaccard similarity
* MinHash vs minimizer similarity
* Runtime behavior

Experiment scripts and recorded results are stored in the
`experiments/` directory.

---

## Real Genome Experiment

The project was tested on real genome sequences including:

* E. coli K-12
* Shigella

For the E. coli vs Shigella comparison using:

```text
k = 15
sketch size = 1000
```

the MinHash implementation produced:

```text
Jaccard similarity: 0.539
Mash distance: 0.02373502549132912
```

### k-mer Size Sweep

|  k | Jaccard Similarity | Mash Distance |
| -: | -----------------: | ------------: |
| 15 |              0.539 |      0.023735 |
| 21 |              0.482 |      0.020479 |
| 31 |              0.407 |      0.017653 |

These experiments show how the estimated similarity changes with the
choice of k-mer size.

---

## Testing

The project currently contains **66 automated tests**.

The tests cover:

* Canonical k-mer generation
* Reverse complements
* Invalid DNA characters
* MinHash sketching
* Deterministic hashing
* Jaccard similarity
* Mash distance
* FASTA parsing
* Minimizer generation
* Minimizer-based sketching
* Minimizer reproducibility
* Minimizer window behavior
* Minimizer similarity trade-offs
* Command-line execution

### Run Tests

```bash
pytest
```

Expected result:

```text
66 passed
```

---

## Continuous Integration

GitHub Actions automatically runs the test suite on pushes and pull
requests.

The workflow uses:

```text
Python 3.11
```

and installs dependencies from:

```text
requirements.txt
```

---

## Reproducibility

The implementation uses deterministic MD5 hashing rather than Python's
built-in `hash()` function.

This ensures that hash-based results are reproducible across different
Python processes.

Experiments are stored as scripts and recorded results in the
`experiments/` directory.

---

## Project Goals

The long-term goal of this project is to build a deeper understanding of
the algorithms and computational techniques used in large-scale genome
comparison.

The project is being developed progressively from basic sequence
processing toward more efficient and research-oriented genome similarity
methods.

Future directions include:

* Improved minimizer strategies
* Weighted minimizers
* Repeat-aware sampling
* More extensive benchmarking
* Accuracy comparison against established tools
* Memory profiling
* Larger genome datasets

---

## License

This project is released under the MIT License.
