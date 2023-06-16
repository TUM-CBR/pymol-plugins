from typing import Dict

CLUSTAL_SPACERS = ['-']

def is_spacer(item : str):
    return item in CLUSTAL_SPACERS

def remove_spacers(sequence : str):

    for spacer in CLUSTAL_SPACERS:
        sequence = sequence.replace(spacer, '')

    return sequence

def parse_alignments(alignments : str) -> Dict[str, str]:

    results = {}

    for line in alignments.splitlines():
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