from typing import Iterator

from Bio.Seq import Seq

class SeqRecord:
    id : str
    _seq : Seq

    def __iter__(self) -> Iterator[str]: ...