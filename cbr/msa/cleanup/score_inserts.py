from Bio.Align import MultipleSeqAlignment
from typing import List, Optional

from ...clustal.msa import is_blank
from .score_gap_divergence import gaps_by_position
from .utils import normalized

def score_insert_line(line : List[bool]) -> float:
    score = 0
    ix = 0
    acc = 0

    while(ix < len(line)):
        if not line[ix]:
            score += 0 if acc < 3 else 2**acc
            acc = 0
        else:
            acc+=1
        
        ix += 1

    return score

def score_inserts(
        sequences : MultipleSeqAlignment,
        continuity_treshold : float,
        gaps_in_position : Optional[List[int]] = None
) -> List[float]:
    
    if gaps_in_position is None:
        gaps_in_position = gaps_by_position(sequences)

    inserts = [
        [
            not is_blank(res) and gaps_in_position[ix] / len(sequences) > continuity_treshold
            for ix,res in enumerate(seq)
        ]
        for seq in sequences
    ]

    scores_raw = [
        score_insert_line(line) for line in inserts
    ]

    return normalized(scores_raw)

