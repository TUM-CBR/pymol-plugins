from Bio.Align import MultipleSeqAlignment
import numpy as np
from numpy.typing import NDArray
from typing import Iterable, List, NamedTuple, Set, Tuple

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
    values : Iterable[bool],
    min_segment_size: int = 1
    ) -> List[int]:
    
    segments: Set[Tuple[int,int]] = set()
    segment_start = -1
    values_len = 0

    for ix,value in enumerate(values):
        
        segment_size = ix - segment_start - 1

        if not value and segment_size >= min_segment_size:

            # segment start is set when a 'false' value is
            # observed. Therefore, it's index should not be
            # included in the segment
            segments.add((segment_start + 1, ix))
            segment_start = ix
        elif not value:
            segment_start = ix

        values_len += 1

    result = list(int(0) for _ in range(0, values_len))

    for (start, end) in segments:
        for ix in range(start, end):
            result[ix] = end - start

    return result