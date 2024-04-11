import numpy as np
from numpy.typing import NDArray
from typing import Optional

from ..cleanup.support import score_contigous
from .score_gap_divergence import gaps_by_position
from .support import ScoreContext

def enumerate_ravines(
    continuity_treshold : float,
    position_to_gaps_mapping : NDArray[np.int64],
    ) -> NDArray[np.int64]:

    n_seqs = len(position_to_gaps_mapping)
    continuity_mask = (position_to_gaps_mapping / n_seqs) > continuity_treshold

    return score_contigous(continuity_mask)

def position_to_ravines(
    sequences : ScoreContext,
    continuity_treshold : float,
    position_to_gaps_mapping : Optional[NDArray[np.int64]] = None
) -> NDArray[np.int64]:
    
    if position_to_gaps_mapping is None:
        position_to_gaps_mapping = gaps_by_position(sequences)

    return enumerate_ravines(continuity_treshold, position_to_gaps_mapping)

def score_residues_in_ravines(
        sequences : ScoreContext,
        continuity_treshold : float,
        position_to_gaps_mapping : Optional[NDArray[np.int64]] = None
) -> NDArray[np.int64]:

    ravines = position_to_ravines(sequences, continuity_treshold, position_to_gaps_mapping)
    assert len(ravines) == sequences.alignment.get_alignment_length()
    result = np.zeros(sequences.vectorized.shape, dtype=np.int64)
    np.putmask(result, sequences.vectorized != '-', ravines)

    return result

    

    