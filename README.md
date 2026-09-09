## Current Implementation
![Tests](https://img.shields.io/badge/tests-65%20passed-brightgreen)
[![Tests](https://github.com/Nitin9775/genome-similarity-from-scratch/actions/workflows/tests.yml/badge.svg)](https://github.com/Nitin9775/genome-similarity-from-scratch/actions/workflows/tests.yml)

The project currently supports:

- Canonical k-mer generation
- FASTA file reading
- Deterministic hashing
- Bottom-k MinHash sketches
- Jaccard similarity estimation
- Mash distance calculation

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