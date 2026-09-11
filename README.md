# Genome Similarity From Scratch

A research-oriented implementation of genome similarity estimation using
canonical k-mers, MinHash sketches, Mash distance, streaming minimizers,
and experimental repeat-aware weighting.

The goal of this project is to understand and implement the core algorithms
behind genome similarity estimation rather than simply relying on existing
bioinformatics tools.

---

## Current Implementation

![Tests](https://img.shields.io/badge/tests-74%20passed-brightgreen)

[![Tests](https://github.com/Nitin9775/genome-similarity-from-scratch/actions/workflows/tests.yml/badge.svg)](https://github.com/Nitin9775/genome-similarity-from-scratch/actions/workflows/tests.yml)

The project currently supports:

- Canonical k-mer generation
- Reverse-complement handling
- FASTA file parsing
- Deterministic MD5 hashing
- Bottom-k MinHash sketches
- Jaccard similarity estimation
- Exact Jaccard calculation for evaluation
- Mash distance calculation
- Streaming minimizer generation
- Deque-based minimizer selection
- Minimizer-based genome sketching
- Experimental complexity-weighted minimizers
- Parameter-sweep experiments
- Runtime benchmarking
- Comparison against external Mash 2.3
- Automated testing
- GitHub Actions continuous integration

---

## Project Structure

```text
genome-similarity-from-scratch/
│
├── data/
│   ├── ecoli.fasta
│   ├── shigella.fasta
│   ├── bacillus_subtilis.fasta
│   ├── bacillus_cereus.fasta
│   ├── test.fasta
│   └── test_genome_2.fasta
│
├── experiments/
│   ├── minhash_benchmark.py
│   ├── minhash_accuracy.py
│   ├── minhash_accuracy_large.py
│   ├── minhash_accuracy_sweep.py
│   ├── kmer_runtime_benchmark.py
│   ├── minimizer_sweep.py
│   ├── weighted_minimizer_benchmark.py
│   ├── weighted_minimizer_accuracy.py
│   ├── exact_vs_minhash.py
│   ├── ecoli_shigella.txt
│   ├── minimizer_sweep.txt
│   ├── mash_benchmark.txt
│   ├── method_comparison.txt
│   ├── weighted_minimizer_results.txt
│   ├── bacillus_comparison.txt
│   ├── summary_results.txt
│   └── experiment_config.txt
│
├── src/
│   ├── __init__.py
│   ├── ani_estimator.py
│   ├── cli.py
│   ├── minimizers.py
│   ├── sketching.py
│   └── weighted_minimizers.py
│
├── tests/
│   ├── test_sketching.py
│   └── test_weighted_minimizers.py
│
├── .github/
│   └── workflows/
│       └── tests.yml
│
├── requirements.txt
├── README.md
└── LICENSE
```

---

## Methods

### Canonical k-mers

A DNA sequence is divided into overlapping k-mers.

For each k-mer, both the forward sequence and its reverse complement are
considered. The lexicographically smaller sequence is selected as the
canonical representation.

This prevents a sequence and its reverse complement from being treated as
different k-mers.

Invalid k-mers containing characters other than `A`, `C`, `G`, and `T`
are skipped.

### MinHash Sketching

The project implements a bottom-k MinHash sketch from scratch.

Each canonical k-mer is converted into a deterministic MD5 hash value.
Only the smallest hash values up to a configurable sketch size are retained.

The resulting sketches are used to estimate Jaccard similarity between
genomes.

### Jaccard Similarity

For two sets of k-mers, Jaccard similarity is defined as:

```text
J(A,B) = |A ∩ B| / |A ∪ B|
```

A value closer to `1` indicates greater similarity, while a value closer
to `0` indicates lower similarity.

The project also implements exact Jaccard calculation, which is used as a
reference when evaluating the MinHash approximation.

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

## Minimizer-based Sketching

The project implements streaming minimizer selection from scratch.

For each window of consecutive k-mers, the minimum canonical k-mer hash
is selected as the minimizer.

A deque-based streaming algorithm is used so that the complete sequence of
k-mer hashes does not need to be materialized before minimizer selection.

For the E. coli genome with `k=15` and `w=10`:

```text
Canonical k-mers:  4,641,638
Minimizers:          816,575
Reduction:            82.41%
```

For Shigella:

```text
Canonical k-mers:  4,828,786
Minimizers:          784,905
Reduction:            83.75%
```

This demonstrates the substantial reduction in sequence representation
achieved by minimizer sampling.

---

## Minimizer Window Experiment

The effect of minimizer window size was evaluated using E. coli and
Shigella.

| Window | E. coli Minimizers | Shigella Minimizers | Jaccard |
| -----: | -----------------: | ------------------: | ------: |
|      5 |          1,494,480 |           1,436,830 |  0.5335 |
|     10 |            816,575 |             784,905 |  0.5278 |
|     20 |            428,206 |             412,133 |  0.5173 |
|     50 |            176,237 |             169,840 |  0.5023 |
|    100 |             89,107 |              85,667 |  0.4943 |

Increasing the window size reduces the number of selected minimizers, but
also changes the observed similarity.

From `w=5` to `w=100`:

```text
Jaccard: 0.5335 → 0.4943
Relative change: -7.36%
```

This demonstrates the trade-off between representation reduction and the
similarity signal retained by minimizer sampling.

---

## MinHash Accuracy Experiments

The implementation was evaluated against exact Jaccard similarity.

For a controlled sequence experiment with sequence length 3,500 and
`k=15`, increasing the sketch size improved the approximation in cases
where the exact Jaccard value was non-zero.

Example:

```text
Exact Jaccard: 0.195312

Sketch size    MinHash       Absolute error
50             0.140000     0.055312
100            0.180000     0.015313
500            0.195312     0.000000
1000           0.195312     0.000000
```

These experiments illustrate the relationship between sketch size and
MinHash approximation accuracy.

---

## Runtime Benchmark

The MinHash implementation was benchmarked on E. coli.

After optimizing the sketch update path using a cached maximum hash:

| Sketch Size | Runtime |
| ----------: | ------: |
|         100 | 14.60 s |
|         500 | 13.11 s |
|        1000 | 13.22 s |
|        5000 | 20.91 s |

The experiments demonstrate how sketch size affects computational cost.

---

## E. coli vs Shigella

Using:

```text
k = 15
MinHash sketch size = 1000
```

the project produced:

```text
Exact Jaccard:      0.539657
Our MinHash:        0.539000
Our Mash distance:  0.0237350
```

Using streaming minimizers with `w=10`:

```text
Minimizer Jaccard:  0.527781
Mash distance:     0.0246495
```

The difference between exact and MinHash Jaccard was small:

```text
Absolute error = 0.000657
```

---

## Comparison Against External Mash

The implementation was compared against Mash 2.3 using the same real
genome pairs.

### E. coli vs Shigella

```text
Our MinHash Mash distance:       0.0237350
Our Minimizer Mash distance:     0.0246495
External Mash 2.3 distance:      0.0192445

External Mash shared hashes:     501/1000
```

The comparison provides an external reference for evaluating the
implementation.

### Bacillus subtilis vs Bacillus cereus

```text
Our MinHash Mash distance:       0.2494420
Our Minimizer Mash distance:     0.2345572
External Mash 2.3 distance:      0.295981

External Mash shared hashes:     1/1000
```

The second genome pair provides an additional evaluation on genomes with
much lower estimated similarity.

---

## Complexity-Weighted Minimizers

An experimental repeat-aware minimizer strategy was implemented using
k-mer sequence complexity.

The goal was to investigate whether simple sequence-complexity weighting
could alter minimizer selection and improve similarity estimation.

For E. coli vs Shigella with `k=15` and `w=10`:

```text
Exact Jaccard:              0.539657

Ordinary minimizer Jaccard: 0.527781
Ordinary error:             0.011876

Weighted minimizer Jaccard: 0.520648
Weighted error:             0.019010
```

The weighting changed minimizer selection:

```text
Ordinary minimizers:  816,575
Weighted minimizers:  844,940
```

However, the weighted method did not improve accuracy or runtime on this
dataset.

This is retained as a negative experimental result rather than being
presented as an improvement. It demonstrates an important aspect of
algorithm development: a plausible weighting strategy does not
necessarily improve downstream similarity estimation.

---

## Command-Line Interface

Two genome FASTA files can be compared directly from the command line.

### MinHash

```bash
python -m src.cli data/ecoli.fasta data/shigella.fasta --k 15 --size 1000
```

Example output:

```text
Jaccard similarity: 0.539
Mash distance: 0.02373502549132912
```

### Minimizers

```bash
python -m src.cli data/ecoli.fasta data/shigella.fasta --minimizers --k 15 --w 10
```

The CLI supports:

```text
--k
--size
--minimizers
--w
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

## Testing

The project currently contains **74 automated tests**.

The tests cover:

* Canonical k-mer generation
* Reverse complements
* Invalid DNA characters
* FASTA parsing
* MinHash sketching
* Deterministic hashing
* Jaccard similarity
* Mash distance
* Exact Jaccard calculation
* Minimizer generation
* Minimizer window behavior
* Minimizer-based sketching
* Minimizer reproducibility
* Weighted minimizer generation
* Weighted minimizer reproducibility
* Command-line execution

Run the test suite with:

```bash
pytest
```

Current result:

```text
74 passed
```

---

## Continuous Integration

GitHub Actions automatically runs the test suite on pushes and pull
requests.

The workflow uses Python 3.11 and installs dependencies from:

```text
requirements.txt
```

---

## Reproducibility

The implementation uses deterministic MD5 hashing rather than Python's
built-in `hash()` function.

This ensures that hash-based results remain reproducible across different
Python processes.

Experimental configurations and results are recorded in the
`experiments/` directory.

The main reproducibility configuration is stored in:

```text
experiments/experiment_config.txt
```

---

## Limitations

This project is an educational and research-oriented implementation rather
than a replacement for production genome comparison tools.

In particular:

* The current MinHash implementation uses MD5 for deterministic hashing.
* The Mash distance is calculated from the project's estimated Jaccard
  similarity.
* Minimizer sampling can introduce approximation differences relative to
  exact k-mer Jaccard similarity.
* The complexity-weighted minimizer strategy tested here did not improve
  accuracy on the evaluated dataset.
* The external Mash comparisons demonstrate differences between this
  implementation and the established Mash 2.3 implementation.

These limitations provide directions for future algorithmic investigation.

---

## Project Goals

The main goal of this project is to build a deeper understanding of the
algorithms and computational techniques used in large-scale genome
comparison.

The implementation was developed progressively from:

```text
DNA sequences
      ↓
Canonical k-mers
      ↓
Deterministic hashing
      ↓
MinHash sketches
      ↓
Jaccard similarity
      ↓
Mash distance
      ↓
Streaming minimizers
      ↓
Parameter experiments
      ↓
External tool comparison
      ↓
Experimental weighted minimizers
```

The project emphasizes understanding, implementation, experimentation,
benchmarking, and reproducibility.

---

## License

This project is released under the MIT License.
