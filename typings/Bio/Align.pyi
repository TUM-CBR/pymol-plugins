import numpy as np
from numpy.typing import NDArray
from typing import Iterator, Sequence

from Bio.SeqRecord import SeqRecord

class Alignment:
    ...

class MultipleSeqAlignment:

    def __iter__(self) -> Iterator[SeqRecord]: ...

    @property
    def alignment(self) -> Alignment: ...

    @property
    def indices(self) -> NDArray[np.int64]: ...

    @property
    def inverse_indices(self) -> Sequence[NDArray[np.int64]]: ...

    def __getitem__(self, index: int) -> SeqRecord: ...