from collections import deque
import hashlib


def minimizers(seq, k, w):
    """Generate minimizers from a DNA sequence.

    A minimizer is the smallest canonical k-mer hash
    within each window of w consecutive k-mers.
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
            16
        )

        while window and window[0][0] <= i - w:
            window.popleft()

        while window and window[-1][1][0] >= hash_value:
            window.pop()

        window.append((i, (hash_value, canonical_kmer)))

        if i >= min(w, len(seq) - k + 1) - 1 and window:
            current = window[0][1]

            if current != last_minimizer:
                yield current[1]
                last_minimizer = current