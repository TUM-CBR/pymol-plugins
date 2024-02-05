from typing import Dict

from ..blosum import blosum62, blosum80, BlosumMatrix

BLOSUM_MATRIXES : Dict[str, BlosumMatrix] = {
    "BLOSUM 62": blosum62,
    "BLOSUM 80": blosum80
}