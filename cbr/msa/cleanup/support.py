from Bio.Align import MultipleSeqAlignment
import numpy as np
from numpy.ma import MaskedArray
from numpy.typing import NDArray
from typing import Any, NamedTuple

class ScoreContext(NamedTuple):
    alignment : MultipleSeqAlignment
    vectorized : NDArray[np.uintc]

    @classmethod
    def from_alignment(cls, alignment: MultipleSeqAlignment) -> 'ScoreContext':
        return ScoreContext(
            alignment=alignment,
            vectorized=np.stack([list(seq) for seq in alignment])
        )

def normalized(values : NDArray[np.float64]) -> NDArray[np.float64]:
    min_v = np.min(values)
    max_v = np.max(values)
    spread = max_v - min_v

    if abs(spread) < 1:
        spread = 1

    return (values - min_v) / spread

def ranked(values : NDArray[np.float64]) -> NDArray[np.int64]:
    sorted_indexes = np.argsort(values)
    result = np.zeros(values.shape, dtype=np.int64)
    result[sorted_indexes] = np.arange(0, len(sorted_indexes))
    return result

def score_contigous(
    values : NDArray[np.bool_],
    min_segment_size: int = 1
    ) -> NDArray[np.int64]:
    """
    This function accepts an array of boolean values. Each index of the array corresponds to
    and index in a multiple sequence alignment and True/False indicates wether a property
    is present or absent at that position.

    This function finds the lenghts of the continous segments which has the value True
    and replaces the Trues with the lenght of the segment and False with 0.

    Parameters
    ----------
    values : NDArray[np.bool_]
        The array of True/False values
    min_segment_size : int
        The minimum  number of contigous values for the replacement to be applied. 0 is used when
        this treshold is not reached.
    """

    result = np.zeros(values.shape, dtype=np.int64)
    result_ma : MaskedArray[np.int64, Any] = np.ma.array(result, mask=np.logical_not(values))

    for slice in np.ma.flatnotmasked_contiguous(result_ma):
        start = slice.start
        stop = slice.stop
        value = stop - start

        result[start:stop] = value
    
    return result