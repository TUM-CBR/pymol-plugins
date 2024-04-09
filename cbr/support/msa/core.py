import itertools
from Bio.Align import MultipleSeqAlignment
import numpy  as np
from numpy.typing import NDArray
import statistics
from typing import Iterable, NamedTuple, Dict, List, Sequence, Tuple

Msa = Dict[str, str]

def from_biopython(alignment: MultipleSeqAlignment) -> Msa:
    return dict(
            (seq.id, str(seq._seq))
            for seq in alignment
        )

def enumerate_pairs(msa: Msa, pos_1: int, pos_2: int) -> Dict[Tuple[str, str], int]:

    result : Dict[Tuple[str, str], int] = {}

    for sequence in msa.values():
        pair = (sequence[pos_1].upper(), sequence[pos_2].upper())

        if pair not in result:
            result[pair] = 0

        result[pair] += 1

    return result


class CoEvolving(NamedTuple):
    """
    A class used to represent a set of positions and residues of a multiple
    sequence alignment which occur in a correlated fashion.

    Attributes
    __________
    alignment : MultipleSequenceAlignment
        The MSA from which this analysis was constructed 
    positions : Tuple[int]
        The positions in the multiple sequence alignment where the correlation is observed
    residues : Tuple[str]
        The residue for which the co-evolution was observed at the positions above
    occurence : float
        The percentage of sequences where the combination of residues and possitions occurs
    occurence_ratio : float
        The ratio of sequences that contain the co-evolving residues vs the sequences that contain
        any of the residues.
    symmetry : float
        How frequently do the same residues occur in the different positions captured by this record.
    """

    alignment: MultipleSeqAlignment
    positions : Tuple[int, ...]
    residues : Tuple[str, ...]
    occurence : float
    occurence_ratio : float
    symmetry : float

GAPS = ['-']

def enumerate_coevolving(
    alignment: MultipleSeqAlignment,
    min_treshold: float = 0.1,
    max_treshold: float = 0.9
) -> Iterable[CoEvolving]:
    """
    A function to identify pairs of residues that show co-evolution in a multiple
    sequence alignment.

    This function provides two parameters to perform an initial filtering of the
    co-evolving space. The "min_treshold" parameter considers the minimum percentage
    of sequences where the residues must be observed in order to be considered a
    co-evolving pair. The "max_treshold" parameter consideres the maximum percentage
    or sequences where the residues can be observed. In other words, exceding this
    treshold would consider the resiude set as 'conserved'.

    Parameters
    __________
    alignment : MultipleSequenceAlignment
        The alignment to be analyzed
    min_treshold : float
        The minimum number of occurences of a set of co-evolving residues to be considered
    max_treshold : float
        The maximum number of occurences of a set of co-evolving residues to be considered
    """

    sequences : NDArray[np.uintc] = np.stack([list(seq) for seq in alignment])
    (n_seqs, msa_len) = sequences.shape

    for i in range(0, msa_len):

        current_column = sequences[:,i]
        current_counts = np.unique(current_column, return_counts=True)

        for (current_res, current_count) in zip(*current_counts):

            if current_res in GAPS or current_count/n_seqs < min_treshold:
                continue

            current_mask = current_column == current_res

            for j in range(i+1, msa_len):

                reference_column = sequences[:,j]
                candidates_counts = np.unique(reference_column[current_mask], return_counts=True)
                symmetry_dict : Dict[Tuple[int,...], float] = {}
                current_results : List[CoEvolving] = []

                for candidate_res, candidate_count in zip(*candidates_counts):

                    occurence = candidate_count / n_seqs
                    if candidate_res in GAPS \
                        or occurence < min_treshold \
                        or occurence > max_treshold:
                        continue

                    candiate_mask = reference_column == candidate_res
                    occurence_ratio = candidate_count / len(reference_column[current_mask | candiate_mask])
                    residues = (current_res, candidate_res)

                    current_results.append(
                        CoEvolving(
                            alignment=alignment,
                            positions=(i,j),
                            residues=residues,
                            occurence=occurence,
                            occurence_ratio=occurence_ratio,
                            symmetry=0
                        )
                    )
                    symmetry_dict[residues] = occurence

                for result in current_results:
                    symmetry = statistics.stdev(
                        value if value is not None else 0
                        for permutation in itertools.permutations(result.positions)
                        for value in [symmetry_dict.get(permutation)]
                    )

                    yield result._replace(symmetry=1 - symmetry)


class CoEvolutionAnalyzer:

    def __init__(self, entries: Sequence[CoEvolving]):

        self.__entries = entries
        pos_index = {}

        self.__alignment = entries[0].alignment

        for entry in entries:

            key = entry.positions