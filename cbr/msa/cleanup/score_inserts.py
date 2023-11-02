from Bio.Align import MultipleSeqAlignment
from typing import List, Optional

from .score_gap_divergence import gaps_by_position

def score_insert_line(line : List[bool]) -> int:
    score = 0
    ix = 0
    acc = 0

    while(ix < len(line)):
        if not line[ix]:
            score += acc
            acc = 0
        elif acc == 0:
            acc = 2
        else:
            acc*=2
        
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
        [gaps_in_position[ix] / len(sequences) < continuity_treshold for ix,_ in enumerate(seq)]
        for seq in sequences
    ]

    scores_raw = [
        score_insert_line(line) for line in inserts
    ]

    max_score = max(scores_raw)

    return [score/max_score for score in scores_raw]

