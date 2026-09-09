from src.sketching import kmerize, MinHashSketch, read_fasta, sketch_fasta, sketch_sequence, sketch_sequence_with_minimizers 
from src.ani_estimator import mash_distance
import pytest
import subprocess
import sys

def test_kmerize():
    sequence = "ATGCAT"
    result = list(kmerize(sequence, 3))

    assert result == ["ATG", "GCA", "GCA", "ATG"]

def test_reverse_complement():
    sequence = "ATGC"
    result = list(kmerize(sequence, 4))

    reverse_complement_result = list(kmerize("GCAT", 4))

    assert result == reverse_complement_result

def test_minhash_identical():
    sketch1 = MinHashSketch(size=10)
    sketch2 = MinHashSketch(size=10)

    for kmer in ["ATG", "GCA", "TTA"]:
        sketch1.add(kmer)
        sketch2.add(kmer)

    assert sketch1.jaccard(sketch2) == 1.0

def test_minhash_different():
    sketch1 = MinHashSketch(size=10)
    sketch2 = MinHashSketch(size=10)

    for kmer in ["ATG", "GCA", "TTA"]:
        sketch1.add(kmer)

    for kmer in ["ATG", "GCA", "CCC"]:
        sketch2.add(kmer)

    assert 0.0 < sketch1.jaccard(sketch2) < 1.0

def test_deterministic_hash():
    sketch1 = MinHashSketch(size=10)
    sketch2 = MinHashSketch(size=10)

    sketch1.add("ATG")
    sketch2.add("ATG")

    assert sketch1.hashes == sketch2.hashes

def test_minhash_from_sequence():
    sketch = MinHashSketch(size=10)
    sketch.from_sequence("ATGCAT", 3)

    assert len(sketch.hashes) == 2

def test_read_fasta():
    sequence = read_fasta("data/test.fasta")

    assert sequence == "ATGCATGCATGCATGCATGC"

def test_sketch_sequence():
    sequence = "ATGCATGCATGCATGCATGC"

    sketch = sketch_sequence(sequence, k=3, size=10)

    assert len(sketch.hashes) > 0
    assert len(sketch.hashes) <= 10

def test_sketch_fasta():
    sketch = sketch_fasta("data/test.fasta", k=3, size=10)

    assert len(sketch.hashes) > 0
    assert len(sketch.hashes) <= 10

def test_genome_similarity():
    sketch1 = sketch_fasta("data/test.fasta", k=3, size=10)
    sketch2 = sketch_fasta("data/test_genome_2.fasta", k=3, size=10)

    similarity = sketch1.jaccard(sketch2)

    assert 0.0 <= similarity <= 1.0
    assert similarity > 0.0
    assert similarity < 1.0

def test_mash_distance():
    distance = mash_distance(1.0, 21)

    assert distance == 0.0

def test_mash_distance_nonzero():
    distance = mash_distance(0.5, 21)

    assert distance > 0.0

def test_mash_distance_ordering():
    low_similarity = mash_distance(0.5, 21)
    high_similarity = mash_distance(0.8, 21)

    assert high_similarity < low_similarity

def test_mash_distance_invalid_jaccard():
    with pytest.raises(ValueError):
        mash_distance(0.0, 21)


def test_mash_distance_invalid_k():
    with pytest.raises(ValueError):
        mash_distance(0.5, 0)

def test_mash_distance_perfect_similarity():
    assert mash_distance(1.0, 15) == 0.0

def test_kmerize_sequence_shorter_than_k():
    result = list(kmerize("ATG", 5))
    assert result == []

def test_kmerize_empty_sequence():
    with pytest.raises(ValueError):
        list(kmerize("", 3))

def test_cli():
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "src.cli",
            "data/test.fasta",
            "data/test_genome_2.fasta",
            "--k",
            "3",
            "--size",
            "10",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert "Jaccard similarity:" in result.stdout
    assert "Mash distance:" in result.stdout

def test_invalid_sketch_size():
    with pytest.raises(ValueError):
        MinHashSketch(size=0)

def test_identical_genomes_have_perfect_similarity():
    sketch1 = sketch_fasta("data/test.fasta", k=3, size=10)
    sketch2 = sketch_fasta("data/test.fasta", k=3, size=10)

    assert sketch1.jaccard(sketch2) == 1.0

def test_empty_sketch_similarity():
    sketch1 = MinHashSketch(size=10)
    sketch2 = MinHashSketch(size=10)

    assert sketch1.jaccard(sketch2) == 0.0

def test_kmerize_invalid_k():
    with pytest.raises(ValueError):
        list(kmerize("ATGC", 0))

def test_kmerize_k_larger_than_sequence():
    result = list(kmerize("ATGC", 10))
    assert result == []

def test_duplicate_kmers_do_not_duplicate_hashes():
    sketch = MinHashSketch(size=10)

    sketch.add("ATG")
    sketch.add("ATG")
    sketch.add("ATG")

    assert len(sketch.hashes) == 1

def test_sketch_never_exceeds_size():
    sketch = MinHashSketch(size=5)

    for kmer in ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT"]:
        sketch.add(kmer)

    assert len(sketch.hashes) <= 5

def test_sketch_is_reproducible():
    sequence = "ATGCATGCATGCATGC"

    sketch1 = sketch_sequence(sequence, k=5, size=5)
    sketch2 = sketch_sequence(sequence, k=5, size=5)

    assert sketch1.hashes == sketch2.hashes

def test_different_sketch_sizes():
    sequence = "ATGCATGCATGCATGCATGC"

    small = sketch_sequence(sequence, k=3, size=3)
    large = sketch_sequence(sequence, k=3, size=10)

    assert len(small.hashes) <= 3
    assert len(large.hashes) <= 10
    assert small.hashes.issubset(large.hashes)

def test_small_sketch_is_lowest_hashes():
    sequence = "ATGCATGCATGCATGCATGC"

    sketch = sketch_sequence(sequence, k=3, size=3)

    assert sketch.hashes == set(sorted(sketch.hashes)[:3])

def test_jaccard_is_symmetric():
    sketch1 = sketch_sequence("ATGCATGCATGC", k=3, size=5)
    sketch2 = sketch_sequence("ATGCATGAAAAA", k=3, size=5)

    assert sketch1.jaccard(sketch2) == sketch2.jaccard(sketch1)

def test_jaccard_with_one_empty_sketch():
    sketch1 = MinHashSketch(size=10)
    sketch2 = sketch_sequence("ATGCATGC", k=3, size=10)

    assert sketch1.jaccard(sketch2) == 0.0

def test_jaccard_disjoint_sketches():
    sketch1 = MinHashSketch(size=10)
    sketch2 = MinHashSketch(size=10)

    sketch1.add("AAA")
    sketch2.add("CCC")

    assert sketch1.jaccard(sketch2) == 0.0

def test_jaccard_identical_sketches():
    sketch1 = sketch_sequence("ATGCATGC", k=3, size=10)
    sketch2 = sketch_sequence("ATGCATGC", k=3, size=10)

    assert sketch1.jaccard(sketch2) == 1.0

def test_larger_sketch_contains_smaller_sketch():
    sequence = "ATGCATGCATGCATGCATGCATGC"

    small = sketch_sequence(sequence, k=3, size=5)
    large = sketch_sequence(sequence, k=3, size=10)

    assert small.hashes.issubset(large.hashes)

def test_kmerize_k_equals_one():
    result = list(kmerize("ATGC", 1))

    assert result == ["A", "A", "C", "C"]

def test_kmerize_lowercase_sequence():
    result = list(kmerize("atgc", 3))

    assert result == ["ATG", "GCA"]

def test_kmerize_invalid_nucleotide():
    result = list(kmerize("ATGX", 3))

    assert result == ["ATG"]

def test_kmerize_normalizes_lowercase():
    uppercase = list(kmerize("ATGC", 3))
    lowercase = list(kmerize("atgc", 3))

def test_kmerize_normalizes_lowercase():
    uppercase = list(kmerize("ATGC", 3))
    lowercase = list(kmerize("atgc", 3))

    assert uppercase == lowercase

def test_sketch_sequence_invalid_k():
    with pytest.raises(ValueError):
        sketch_sequence("ATGCATGC", k=0, size=10)

def test_sketch_fasta_invalid_k():
    with pytest.raises(ValueError):
        sketch_fasta("data/test.fasta", k=0, size=10)

def test_sketch_sequence_invalid_size():
    with pytest.raises(ValueError):
        sketch_sequence("ATGCATGC", k=3, size=0)

def test_mash_distance_ordering():
    low_similarity = mash_distance(0.5, 21)
    high_similarity = mash_distance(0.8, 21)

    assert high_similarity < low_similarity

def test_mash_distance_high_similarity():
    distance = mash_distance(0.99, 21)

    assert distance > 0
    assert distance < 0.01

def test_mash_distance_zero_similarity_is_invalid():
    with pytest.raises(ValueError):
        mash_distance(0.0, 21)

def test_mash_distance_invalid_k_negative():
    with pytest.raises(ValueError):
        mash_distance(0.5, -1)

def test_mash_distance_one_similarity():
    assert mash_distance(1.0, 21) == 0.0

def test_jaccard_is_between_zero_and_one():
    sketch1 = sketch_sequence("ATGCATGCATGC", k=3, size=5)
    sketch2 = sketch_sequence("CCCCCCCCCCCC", k=3, size=5)

    similarity = sketch1.jaccard(sketch2)

    assert 0.0 <= similarity <= 1.0

def test_sketch_sequence_empty():
    with pytest.raises(ValueError):
        sketch_sequence("", k=3, size=10)

def test_sketch_fasta_empty():
    with pytest.raises(ValueError):
        sketch_fasta("data/empty.fasta", k=3, size=10)

def test_sketch_fasta_header_only():
    with pytest.raises(ValueError):
        sketch_fasta("data/header_only.fasta", k=3, size=10)

def test_read_fasta_multiple_lines():
    sequence = read_fasta("data/test.fasta")

    assert sequence == "ATGCATGCATGCATGCATGC"

def test_fasta_sequence_can_be_sketches():
    sequence = read_fasta("data/test.fasta")
    sketch = sketch_sequence(sequence, k=3, size=10)

    assert len(sketch.hashes) > 0
    assert len(sketch.hashes) <= 10

from src.minimizers import minimizers


def test_minimizers_basic():
    result = list(minimizers("ATGCATGCATGC", 3, 4))

    assert result == ["ATG"]


def test_minimizers_reverse_complement():
    seq1 = "ATGCATGCATGC"
    seq2 = "GCATGCATGCAT"

    assert list(minimizers(seq1, 3, 4)) == list(
        minimizers(seq2, 3, 4)
    )


def test_minimizers_invalid_k():
    try:
        list(minimizers("ATGC", 0, 2))
        assert False
    except ValueError:
        assert True


def test_minimizers_invalid_window():
    try:
        list(minimizers("ATGC", 3, 0))
        assert False
    except ValueError:
        assert True

def test_minimizers_reproducible():
    seq = "ATGCATGCATGCATGC"

    result1 = list(minimizers(seq, 3, 4))
    result2 = list(minimizers(seq, 3, 4))

    assert result1 == result2

def test_minimizers_larger_sequence():
    seq = "ATGCGTACGTTAGCGATCGATCGTACG"

    result = list(minimizers(seq, 5, 4))

    assert result == [
        "GCGTA",
        "AACGT",
        "CGTTA",
        "AGCGA",
        "CGATC",
        "ACGAT",
    ]

def test_sketch_sequence_with_minimizers():
    sequence = "ATGCGTACGTTAGCGATCGATCGTACG"

    sketch = sketch_sequence_with_minimizers(sequence, k=5, w=4)

    assert len(sketch.hashes) == 6