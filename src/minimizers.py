from collections import deque


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

    def canonical(kmer):
        complement = str.maketrans("ACGT", "TGCA")
        reverse = kmer.translate(complement)[::-1]
        return min(kmer, reverse)

    hashes = []

    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]

        if any(base not in "ACGT" for base in kmer):
            hashes.append(None)
            continue

        canonical_kmer = canonical(kmer)
        hash_value = hash(canonical_kmer)
        hashes.append((hash_value, canonical_kmer))

    w = min(w, len(hashes))

    window = deque()
    last_minimizer = None

    for i, value in enumerate(hashes):
        while window and window[0][0] <= i - w:
            window.popleft()

        if value is not None:
            while window and window[-1][1][0] >= value[0]:
                window.pop()

            window.append((i, value))

        if i >= w - 1 and window:
            current = window[0][1]

            if current != last_minimizer:
                yield current[1]
                last_minimizer = current