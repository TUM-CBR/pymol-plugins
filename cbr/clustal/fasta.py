from io import TextIOBase
import os
from typing import Iterable, Tuple, TypeVar

from ..core.WrapIO import WrapIO

FastaInput = TypeVar('FastaInput', str, Iterable[str], TextIOBase)

allowed_characters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', 'B']

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

def parse_fasta(in_fasta : FastaInput) -> 'Iterable[Tuple[str, str] | Exception]':

    result = []
    def yield_(k, v = None):

        if v is None:
            result.append(k)
        else:
            result.append((k,v))

    with get_fasta_input_arg(in_fasta) as in_fasta_stream:
        current_header : 'None | str' = None
        current_seq : 'None | str' = None

        for line in in_fasta_stream.stream.readlines():
            line = line.strip()

            if line.startswith(">"):
                if current_header and current_seq:
                    yield_(current_header, current_seq)

                current_header = line[1:]
                current_seq = ''
                continue
            
            line = line.replace(" ", "")
            if current_seq is not None and \
                all([c in allowed_characters for c in line.upper()]):
                current_seq += line
            elif current_header:
                current_header = None
                current_seq = None
                yield_(ValueError('Could not parse line "%s"' % line))

        if current_header and current_seq:
            yield_(current_header, current_seq)

    return result