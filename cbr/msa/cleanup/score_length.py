import numpy as np
from numpy.typing import NDArray

from .support import ScoreContext

def score_length(sequences : ScoreContext) -> NDArray[np.int64]:

    msa = sequences.vectorized
    result = np.zeros(msa.shape, dtype=np.int64)
    np.putmask(result, msa != '-', 1)
    return np.sum(result, axis=1)