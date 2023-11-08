from Bio.Align import MultipleSeqAlignment, SeqRecord
from typing import List

from ...clustal import msa
from .utils import normalized

def gaps_by_position(alignment : MultipleSeqAlignment) -> List[int]:
    positions = alignment.get_alignment_length()
    result = [0 for _ in range(0, positions)]

    for sequence in alignment:
        for (i, c) in enumerate(sequence):
            if msa.is_blank(c):
                result[i] += 1

    return result

GapDivergenceScore = List[float]

def score_by_gap_divergence(alignment : MultipleSeqAlignment) -> GapDivergenceScore:
    score_by_position = gaps_by_position(alignment)

    def score_sequence(seq : SeqRecord) -> float:
        return sum(
            score_by_position[i]
            for (i, c) in enumerate(seq) if not msa.is_blank(c)
        )

    abs_scores = [score_sequence(s) for s in alignment]
    return normalized(abs_scores)