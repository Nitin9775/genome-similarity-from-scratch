from src.weighted_minimizers import weighted_minimizers


def test_weighted_minimizers_are_reproducible():
    sequence = "AAAAAAAAAATGCGTACGTAGCTAGCTAGCTAGCTA"

    result1 = list(weighted_minimizers(sequence, k=5, w=5))
    result2 = list(weighted_minimizers(sequence, k=5, w=5))

    assert result1 == result2


def test_weighted_minimizers_return_kmers_of_correct_length():
    sequence = "AAAAAAAAAATGCGTACGTAGCTAGCTAGCTAGCTA"

    result = list(weighted_minimizers(sequence, k=5, w=5))

    assert result
    assert all(len(kmer) == 5 for kmer in result)