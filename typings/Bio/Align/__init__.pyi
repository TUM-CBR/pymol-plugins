from typing import Iterator

from Bio.SeqRecord import SeqRecord

class MultipleSeqAlignment:

    def __iter__(self) -> Iterator[SeqRecord]: ...