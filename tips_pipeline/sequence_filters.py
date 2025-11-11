from typing import List



def extract_confident_subsequences(
    sequence: str,
    aa_scores: List[float],
    min_score: float = 0.75,
    min_length: int = 8,
) -> List[str]:
    """
    Given:
      - sequence: peptide string (AA letters)
      - aa_scores: per-residue confidence scores, same length as sequence
    This function:
      1. Finds contiguous regions where every residue score > min_score.
      2. Keeps only regions length >= min_length.
      3. For each region, also finds the sub-window (>=min_length) with the
         highest average score and returns that as well.

    Returns:
      list of high-quality subsequences (strings).
    """
    assert len(sequence) == len(aa_scores), \
        "sequence and aa_scores length mismatch"

    result = []
    n = len(sequence)
    i = 0
    while i < n:
        if aa_scores[i] > min_score:
            start = i
            while i < n and aa_scores[i] > min_score:
                i += 1
            end = i

            if end - start >= min_length:
                subseq = sequence[start:end]
                subseq_scores = aa_scores[start:end]
                result.append(subseq)

                # find best-scoring window (length>=min_length)
                best_avg = -1.0
                best_sub = ""
                L = len(subseq)
                for j in range(L - min_length + 1):
                    for k in range(j + min_length, L + 1):
                        window_scores = subseq_scores[j:k]
                        avg_score = sum(window_scores) / len(window_scores)
                        if avg_score > best_avg:
                            best_avg = avg_score
                            best_sub = subseq[j:k]
                result.append(best_sub)
        else:
            i += 1
    return result
