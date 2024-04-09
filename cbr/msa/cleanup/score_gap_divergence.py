import numpy as np
from numpy.typing import NDArray

from .support import ScoreContext

def gaps_by_position(context : ScoreContext) -> NDArray[np.int64]:

    msa = context.vectorized
    result = np.zeros(msa.shape, dtype=np.int64)
    gap_mask = msa == '-'
    np.putmask(result, gap_mask, 1)
    np.putmask(result, gap_mask == False, 0)
    return np.sum(result, axis=0)

def score_by_gap_divergence(alignment : ScoreContext) -> NDArray[np.int64]:
    score_by_position = gaps_by_position(alignment)
    result = np.zeros(alignment.vectorized.shape, dtype=np.int64)
    np.putmask(result, alignment.vectorized == '-', score_by_position)
    return result