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
    """ This function scores a sequence by adding a penalty for each insert
    contained in it. The penalty penalty for each insert is N^2 where N is the
    length of the insert. This is because each residue that is part of the insert
    gets a penalty of value N.

    Parameters:
    position_ravine_mapping: An array indicating what positions have a ravine. This
                             should be the output of calling :func:`position_to_ravines`.
    sequence: The sequence of residues to score.

    Returns:
    An array that matches the indexes of the given sequence and in each position indicates
    the penalty incurred by each of the residues.
    """

    return score_contigous(
        (
            position_ravine_mapping[ix] > 0 and not is_blank(resi)
            for ix,resi in enumerate(sequence)
        ),
        # We require at lest two continous residues to consider
        # something an insert. This can be made a parameter
        # in the future.
        min_segment_size=2
    )

def score_inserts(
        sequences : MultipleSeqAlignment,
        continuity_treshold : float,
        gap_ravines : Optional[List[int]] = None
) -> List[List[int]]:
    """Score the sequences of a Multiple Sequence Alignment using :func:`score_insert_line`
    for each sequence and using the provided alignment to identify the ravines using the
    :func:`position_to_ravines` function.
    """
    
    if gap_ravines is None:
        gap_ravines = position_to_ravines(sequences, continuity_treshold)

    return [
        list(score_insert_line(gap_ravines, seq))
        for seq in sequences
    ]

