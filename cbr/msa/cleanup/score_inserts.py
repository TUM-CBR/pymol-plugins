from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from typing import Iterable, List, Optional

from ...clustal.msa import is_blank
from ..cleanup.support import score_contigous
from .score_gap_ravines import position_to_ravines

def score_insert_line(
    position_ravine_mapping : List[int],
    sequence : SeqRecord
    ) -> Iterable[int]:

    return score_contigous(
        position_ravine_mapping[ix] > 0 and not is_blank(resi)
        for ix,resi in enumerate(sequence)
    )

def score_inserts(
        sequences : MultipleSeqAlignment,
        continuity_treshold : float,
        gap_ravines : Optional[List[int]] = None
) -> List[List[int]]:
    
    if gap_ravines is None:
        gap_ravines = position_to_ravines(sequences, continuity_treshold)

    return [
        list(score_insert_line(gap_ravines, seq))
        for seq in sequences
    ]

