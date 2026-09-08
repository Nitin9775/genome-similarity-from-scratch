import math


def mash_distance(jaccard, k):
    """Calculate Mash distance from Jaccard similarity."""

    if not 0 < jaccard <= 1:
        raise ValueError("Jaccard similarity must be between 0 and 1.")

    if k <= 0:
        raise ValueError("k must be greater than 0.")

    return -(1 / k) * math.log((2 * jaccard) / (1 + jaccard))