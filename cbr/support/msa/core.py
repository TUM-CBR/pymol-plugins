from typing import Dict, Tuple

Msa = Dict[str, str]

def enumerate_pairs(msa: Msa, pos_1: int, pos_2: int) -> Dict[Tuple[str, str], int]:

    result : Dict[Tuple[str, str], int] = {}

    for sequence in msa.values():
        pair = (sequence[pos_1].upper(), sequence[pos_2].upper())

        if pair not in result:
            result[pair] = 0

        result[pair] += 1

    return result
