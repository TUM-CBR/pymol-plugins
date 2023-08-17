from io import TextIOBase
import os
from typing import Dict, Iterable, TypeVar

from ..core.WrapIO import WrapIO

MsaInput = TypeVar('MsaInput', str, TextIOBase)

CLUSTAL_SPACERS = ['-']

def is_spacer(item : str):
    return item in CLUSTAL_SPACERS

def remove_spacers(sequence : str):

    for spacer in CLUSTAL_SPACERS:
        sequence = sequence.replace(spacer, '')

    return sequence

def get_alignments_input(arg : MsaInput) -> WrapIO:

    if isinstance(arg, TextIOBase):
        return WrapIO(stream=arg)

    elif isinstance(arg, str) and os.path.exists(arg):
        return WrapIO(open_stream=arg)
    else:
        return WrapIO(init=arg)

def parse_alignments(in_alignments : MsaInput) -> Dict[str, str]:

    results = {}

    with get_alignments_input(in_alignments) as alignments:

        for line in alignments.stream.readlines():
            entry = line.split()[0:2]

            if len(entry) == 2 and \
                entry[0] != "CLUSTAL" and \
                line[0] not in [' ', '\t']:
                [name, sequence] = entry

                if name in results:
                    results[name] += sequence
                else:
                    results[name] = sequence

    return results

MSA_BLANKS = ["-"]

def is_blank(aa : str) -> bool:
    return aa in MSA_BLANKS

def clean_msa_blanks(msa_sequence : str) -> str:
    for blank in MSA_BLANKS:
        msa_sequence = msa_sequence.replace(blank, "")
    return msa_sequence

def get_relative_positions(msa : Dict[str, str], join_msa : Dict[str, str]) -> Iterable[int]:
    """
    This is a function often used to compute positions in an MSA file relative to an arbitrary
    sequence. The idea is that two inputs are provided, an MSA with multiple sequences and an
    MSA of one of the sequences of the original msa and the sequence one wishes to compute the
    relative positions with. It will return a list where each of the list's indexes can
    be regarded as a position in the original msa and the value of each index can be regarded
    as the position that would correspond to the target sequence. If the value is None,
    it means the position in the MSA file doesn't exactly match any position of the target
    sequence.
    """
    link = next((x for x in join_msa.keys() if x in msa.keys()))
    target = next((x for x in join_msa.keys() if x != link))
    link_seq_all = msa[link]
    link_seq_join = join_msa[link]
    target_seq_join = join_msa[target]

    for i in range(len(link_seq_all)):
        link_fragment_gaps = link_seq_all[0:i]
        link_fragment = clean_msa_blanks(link_fragment_gaps)
        canary = True

        for j in range(len(link_seq_join)):
            join_fragment = link_seq_join[0:j]
            if link_fragment == clean_msa_blanks(join_fragment):
                target_fragment = target_seq_join[0:j]
                ix = len(clean_msa_blanks(target_fragment)) - 1
                yield ix
                canary = False
                break

        if canary:
            raise Exception("The position %i of the msa could not be matched" % i)