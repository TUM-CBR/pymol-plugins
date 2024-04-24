from Bio.Align import MultipleSeqAlignment
from typing import Iterable, Literal, TextIO, Union

Handle = Union[str, TextIO]

AlignFormat = Union[Literal['fasta'], Literal['clustal']]

def write(alignments: Union[MultipleSeqAlignment, Iterable[MultipleSeqAlignment]], handle: Handle, format: AlignFormat) -> int: ...