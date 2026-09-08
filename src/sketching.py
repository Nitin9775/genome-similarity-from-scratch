import hashlib

def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


def kmerize(seq, k):
    """Generate canonical k-mers of length k."""
    
    if not seq:
        raise ValueError("Sequence must not be empty.")

    if k <= 0:
        raise ValueError("k must be greater than 0.")

    seq = seq.upper()

    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]

        if any(base not in "ACGT" for base in kmer):
            continue

        reverse = reverse_complement(kmer)
        yield min(kmer, reverse)

def read_fasta(filepath):
    sequence = []

    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                continue

            sequence.append(line)

    return "".join(sequence)

def sketch_sequence(sequence, k=15, size=1000):
    sketch = MinHashSketch(size=size)
    sketch.from_sequence(sequence, k)
    return sketch

def sketch_fasta(filepath, k=15, size=1000):
    sequence = read_fasta(filepath)
    return sketch_sequence(sequence, k=k, size=size)

class MinHashSketch:
    def __init__(self, size=1000):
        if size <= 0:
            raise ValueError("Sketch size must be greater than 0.")

        self.size = size
        self.hashes = set()

    def add(self, kmer):
        hash_value = int(hashlib.md5(kmer.encode()).hexdigest(), 16)

        if len(self.hashes) < self.size:
            self.hashes.add(hash_value)
            return

        largest_hash = max(self.hashes)

        if hash_value < largest_hash:
            self.hashes.remove(largest_hash)
            self.hashes.add(hash_value)

    def from_sequence(self, seq, k):
        for kmer in kmerize(seq, k):
            self.add(kmer)

    def jaccard(self, other):
        if not self.hashes or not other.hashes:
            return 0.0

        intersection = self.hashes & other.hashes
        union = self.hashes | other.hashes

        return len(intersection) / len(union)