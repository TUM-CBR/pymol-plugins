from Bio.Align import MultipleSeqAlignment
from typing import Iterable, List, Optional

from ...clustal.msa import is_blank
from ..cleanup.support import score_contigous
from .score_gap_divergence import gaps_by_position

def enumerate_ravines(
    continuity_treshold : float,
    position_to_gaps_mapping : List[int],
    sequences: MultipleSeqAlignment
    ) -> Iterable[int]:

    return score_contigous(
        gap_count/len(sequences) > continuity_treshold
        for gap_count in position_to_gaps_mapping
    )

def position_to_ravines(
    sequences : MultipleSeqAlignment,
    continuity_treshold : float,
    position_to_gaps_mapping : Optional[List[int]] = None
):
    if position_to_gaps_mapping is None:
        position_to_gaps_mapping = gaps_by_position(sequences)

    return list(enumerate_ravines(continuity_treshold, position_to_gaps_mapping, sequences))

def score_residues_in_ravines(
        sequences : MultipleSeqAlignment,
        continuity_treshold : float,
        position_to_gaps_mapping : Optional[List[int]] = None
) -> List[List[int]]:

    ravines = position_to_ravines(sequences, continuity_treshold, position_to_gaps_mapping)
    assert len(ravines) == sequences.get_alignment_length()
    return [
        [0 if is_blank(resi) else ravines[i] for i, resi in enumerate(seq)]
        for seq in sequences
    ]

    

    