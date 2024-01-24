from io import TextIOBase
import os
from typing import Dict, Iterable, List, NamedTuple, Tuple, TypeVar, Union

from ..core.WrapIO import WrapIO

FastaInput = TypeVar('FastaInput', str, Iterable[str], TextIOBase)

allowed_characters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', 'B', "-"]

def get_fasta_input_arg(in_fasta : FastaInput) -> WrapIO:

    if isinstance(in_fasta, TextIOBase):
        return WrapIO(stream=in_fasta)
    elif isinstance(in_fasta, str) and os.path.exists(in_fasta):
        return WrapIO(open_stream=in_fasta)

    def init(stream : TextIOBase):

        if isinstance(in_fasta, str):
            stream.write(in_fasta)
        elif isinstance(in_fasta, Iterable):
            for in_value in in_fasta:
                stream.write(in_value)
        else:
            raise ValueError("The argument is not a valid fasta input type")

    return WrapIO(init=init)

def parse_fasta_iter(input_any : FastaInput):
    return parse_fasta(input_any)

def parse_fasta_stream(input: TextIOBase):
    return parse_fasta(input)

class FastaSequenceWithMetadata(NamedTuple):
    id : str
    seq : str
    lines : List[int]

    def as_tuple(self):
        return (self.id, self.seq)

class FastaSequences(NamedTuple):
    sequences : List[FastaSequenceWithMetadata]
    exceptions : List[Exception]

    def as_tuples(self) -> Iterable[Union[Tuple[str, str], Exception]]:
        for seq in self.sequences:
            yield seq.as_tuple()

        for exn in self.exceptions:
            yield exn

    def as_dict(self, throw_on_errors : bool = True) -> Dict[str, str]:

        if throw_on_errors:
            for exn in self.exceptions:
                raise exn

        return dict(
            entry.as_tuple() for entry in self.sequences
        )

def parse_fast_meta(in_fasta : FastaInput) -> FastaSequences:

    sequences : List[FastaSequenceWithMetadata] = []
    errors : List[Exception] = []

    with get_fasta_input_arg(in_fasta) as in_fasta_stream:
        current_header : 'None | str' = None
        current_seq : 'None | str' = None
        current_lines : List[int] = []

        for line in in_fasta_stream.stream.readlines():
            line = line.strip()

            if line.startswith(">"):
                if current_header and current_seq:
                    sequences.append(
                        FastaSequenceWithMetadata(
                            id = current_header,
                            seq = current_seq,
                            lines = current_lines
                        )
                    )

                current_header = line[1:]
                current_seq = ''
                current_lines = []
                continue
            
            line = line.replace(" ", "")
            if current_seq is not None and \
                all([c in allowed_characters for c in line.upper()]):
                current_seq += line

                if len(line) > 0:
                    current_lines.append(len(line))
            elif current_header:
                current_header = None
                current_seq = None
                errors.append(ValueError('Could not parse line "%s"' % line))

        if current_header and current_seq:
            sequences.append(
                FastaSequenceWithMetadata(
                    id = current_header,
                    seq = current_seq,
                    lines = current_lines
                )
            )

    return FastaSequences(
        sequences = sequences,
        exceptions = errors
    )

def parse_fasta(in_fasta : FastaInput) -> 'Iterable[Tuple[str, str] | Exception]':

    return parse_fast_meta(in_fasta).as_tuples()