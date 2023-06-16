from typing import Dict, Iterable, Generator, Tuple

allowed_characters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', 'B']

def parse_fasta_iter(input : str) -> 'Iterable[Tuple[str, str] | Exception]':

    current_header : 'None | str' = None
    current_seq : 'None | str' = None

    for line in input.splitlines():
        line = line.strip()

        if line.startswith(">"):
            if current_header and current_seq:
                yield (current_header, current_seq)

            current_header = line[1:]
            current_seq = ''
        elif current_seq is not None and \
            all([c in allowed_characters for c in line.upper()]):
            current_seq += line
        elif current_header:
            current_header = None
            current_seq = None
            yield ValueError('Could not parse line "%s"' % line)

    if current_header and current_seq:
        yield (current_header, current_seq)