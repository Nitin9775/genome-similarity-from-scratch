from collections import deque
import hashlib


def weighted_minimizers(seq, k, w):
    """Generate minimizers using a simple repeat-aware weight.

    Lower-complexity k-mers receive a higher score, making them
    less likely to be selected than high-complexity k-mers.
    """

    if not seq:
        raise ValueError("Sequence must not be empty.")

    if k <= 0:
        raise ValueError("k must be greater than 0.")

    if w <= 0:
        raise ValueError("Window size must be greater than 0.")

    seq = seq.upper()

    if len(seq) < k:
        return

    complement = str.maketrans("ACGT", "TGCA")

    def canonical(kmer):
        reverse = kmer.translate(complement)[::-1]
        return min(kmer, reverse)

    def score(kmer):
        """Return a repeat-aware score for a canonical k-mer."""
        unique_bases = len(set(kmer))
        return unique_bases

    window = deque()
    last_minimizer = None

    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]

        if any(base not in "ACGT" for base in kmer):
            window.clear()
            last_minimizer = None
            continue

        canonical_kmer = canonical(kmer)
        hash_value = int(
            hashlib.md5(canonical_kmer.encode()).hexdigest(),
            16,
        )

        weighted_hash = hash_value // score(canonical_kmer)

        weighted_value = (
            weighted_hash,
            hash_value,
            canonical_kmer,
        )

        while window and window[0][0] <= i - w:
            window.popleft()

        while window and window[-1][1] >= weighted_value:
            window.pop()

        window.append((i, weighted_value))

        if i >= min(w, len(seq) - k + 1) - 1 and window:
            current = window[0][1]

            if current != last_minimizer:
                yield current[2]
                last_minimizer = current