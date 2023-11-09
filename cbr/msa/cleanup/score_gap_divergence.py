from Bio.Align import MultipleSeqAlignment
from typing import List

from ...clustal import msa

def gaps_by_position(alignment : MultipleSeqAlignment) -> List[int]:
    positions = alignment.get_alignment_length()
    result = [0 for _ in range(0, positions)]

    for sequence in alignment:
        for (i, c) in enumerate(sequence):
            if msa.is_blank(c):
                result[i] += 1

    return result

def score_by_gap_divergence(alignment : MultipleSeqAlignment) -> List[List[int]]:
    score_by_position = gaps_by_position(alignment)

    return [
        [0 if msa.is_blank(resi) else score_by_position[i] for i,resi in enumerate(seq)]
        for seq in alignment
    ]