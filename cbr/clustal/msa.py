import os
from typing import Dict, TextIO, TypeVar

from ..core.WrapIO import WrapIO

MsaInput = TypeVar('MsaInput', str, TextIO)

CLUSTAL_SPACERS = ['-']

def is_spacer(item : str):
    return item in CLUSTAL_SPACERS

def remove_spacers(sequence : str):

    for spacer in CLUSTAL_SPACERS:
        sequence = sequence.replace(spacer, '')

    return sequence

def get_alignments_input(arg : MsaInput) -> WrapIO:

    if isinstance(arg, TextIO):
        return WrapIO(stream=arg)

    elif isinstance(arg, str) and os.path.exists(arg):
        return WrapIO(open_stream=arg)
    else:
        return WrapIO()

def parse_alignments(in_alignments : MsaInput) -> Dict[str, str]:

    results = {}

    with get_alignments_input(in_alignments) as alignments:

        for line in alignments.stream.readlines():
            entry = line.split()

            if len(entry) == 2 and \
                entry[0] != "CLUSTAL" and \
                line[0] not in [' ', '\t']:
                [name, sequence] = entry

                if name in results:
                    results[name] += sequence
                else:
                    results[name] = sequence

    return results