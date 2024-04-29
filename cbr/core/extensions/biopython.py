from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

import numpy as np

def consensus(msa: MultipleSeqAlignment) -> Seq:

    freqs = msa.alignment.frequencies
    res_keys = list(k for k in freqs.keys() if k != "-")

    stacked_freqs = np.vstack([freqs[k] for k in res_keys]).T
    i_max = np.argmax(stacked_freqs, axis=1)

    return Seq(
        "".join(res_keys[i] for i in i_max)
    )
