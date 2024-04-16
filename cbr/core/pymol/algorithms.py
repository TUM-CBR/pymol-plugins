from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from io import StringIO
import numpy as np
from numpy import linalg
from numpy.typing import NDArray
from pymol import cmd
from typing import Any, Dict, cast, Generic, List, NamedTuple, Optional, Tuple, TypeVar

from .structure import StructureSelection

TValue = TypeVar('TValue', np.float16, np.float32, np.float64)

class StructureVector(Generic[TValue]):
    selections: List[StructureSelection]
    indexes: List[NDArray[np.int64]]
    array: NDArray[TValue]

    def __init__(
        self,
        selections: List[StructureSelection],
        indexes: List[NDArray[np.int64]],
        array: NDArray[TValue]
    ):
        self.selections = selections
        self.indexes = indexes
        self.array = array

    def resv_to_position(self, ix: int) -> Dict[int, NDArray[np.float64]]:

        array = self.array[ix] if len(self.array.shape) == 4 else self.array

        return {
            resv: array[i,0]
            for i,resv in enumerate(self.indexes[ix])
        }


ERROR_ATOMS = ["CA"]
ERROR_ATOMS_N = len(ERROR_ATOMS)

def cx_coords(
    s: StructureSelection,
    state: int = 1
) -> StructureVector[np.float64]:
    
    resv_maping:List[int] = []
    ca_coords: List[List[List[float]]] = []

    def apply(resv: int, x: float, y: float, z: float):

        if len(resv_maping) == 0 or resv_maping[-1] != resv:

            while len(ca_coords) > 0 and len(ca_coords[-1]) < ERROR_ATOMS_N:
                ca_coords[-1].append([0,0,0])

            resv_maping.append(resv)
            ca_coords.append([])

        # If we have more than one CA atom, it means
        # there are alternative conformations
        if len(ca_coords[-1]) < ERROR_ATOMS_N:
            ca_coords[-1].append([x, y, z])

    atom_selection = " and ".join(f"name {a}" for a in ERROR_ATOMS)

    cmd.iterate_state(
        state,
        f"{s.base_query} and ({atom_selection})",
        'apply(resv, x, y, z)',
        space={'apply': apply}
    )

    return StructureVector(
        selections=[s],
        indexes=[np.array(resv_maping)],
        array=np.array(ca_coords)
    )

def cx_distance_matrix(
    s1: StructureSelection,
    s2: StructureSelection,
    state1: int = 1,
    state2: int = 1
) -> StructureVector[np.float64]:
    c1 = cx_coords(s1, state1)
    c2 = cx_coords(s2, state2)
    a1 = c1.array
    a2 = c2.array

    deltas = np.array([a2 - a1[i] for i in range(len(a1))])

    norms = np.array([
        [linalg.norm(deltas[i, j]) for j in range(len(a2))]
        for i in range(len(a1))
    ])

    return StructureVector(
        selections=c1.selections + c2.selections,
        indexes=c1.indexes + c2.indexes,
        array=norms
    )

class Interval(NamedTuple):
    bottom: Tuple[int, int]
    top: Tuple[int, int]

    def __valid(self) -> bool:
        (a_low, b_low) = self.bottom
        (a_high, b_high) = self.top

        return a_low <= a_high and b_low <= b_high

    def split(self, index: Tuple[int, int]) -> List['Interval']:
        (ia, ib) = index
        new_items = [
            Interval(bottom=self.bottom, top=(ia-1, ib-1)),
            Interval(bottom=(ia+1, ib+1), top=self.top)
        ]

        return [
            item
            for item in new_items if item.__valid()
        ]

    def contains(self, index: Tuple[int, int]) -> bool:
        (ia, ib) = index
        (a_low, b_low) = self.bottom
        (a_high, b_high) = self.top

        return \
            ia >= a_low and ia <= a_high \
            and ib >= b_low and ib <= b_high

def cx_seq_align(
    distances: StructureVector[np.float64]
):
    d_matrix = distances.array
    ix_iter = np.nditer(d_matrix, flags=['multi_index'])
    indices: List[Tuple[int, int]] = cast(Any, list(ix_iter.multi_index for _ in ix_iter))

    # Sort the indices by the distance in the d_matrix
    # this means that we can iterate the indexes in from closest
    # to furthest
    indices.sort(key=lambda i: d_matrix[i])

    intervals = [Interval(bottom=(0,0), top=(len(d_matrix) - 1, len(d_matrix[0] - 1)))]
    mappings: List[Tuple[int, int]] = []

    for index in indices:

        if len(intervals) == 0:
            break

        # Give the current list of allowed intervals, we check
        # to see if the given index falls into one of them
        candidate = next(
            (
                ix
                for ix,interval in enumerate(intervals)
                    if interval.contains(index)
            ),
            None
        )

        if candidate is None:
            continue

        next_interval = intervals.pop(candidate)
        mappings.append(index)

        for new_interval in next_interval.split(index):
            intervals.append(new_interval)

    mappings.sort(key=lambda i: i[0])

    return mappings

def index_by_resv(index: NDArray[np.int64]):
    return {
        resv: i
        for i,resv in enumerate(index)
    }

class SpatialAlignment(NamedTuple):
    distance_matrix: StructureVector[np.float64]
    alignment: MultipleSeqAlignment
    resv_maps: List[Dict[int, int]]
    resv_to_seq_maps: List[Dict[int, int]]

    def distance_by_resv(self, resv1: int, resv2: int) -> Optional[np.float64]:
        i1 = self.resv_to_seq_maps[0].get(resv1)
        i2 = self.resv_to_seq_maps[1].get(resv2)

        if i1 is None or i2 is None:
            return None

        return self.distance_matrix.array[i1, i2]

def align_by_rmsd(
    s1: StructureSelection,
    s2: StructureSelection,
    state1: int = 1,
    state2: int = 1
) -> SpatialAlignment:
    
    distances = cx_distance_matrix(s1, s2, state1, state2)
    spatial_alignment = cx_seq_align(distances)
    seq1 = s1.get_sequence()
    resvs1 = list(seq1.keys())
    resvs1.sort()
    n_resvs1 = len(resvs1)

    seq2 = s2.get_sequence()
    resvs2 = list(seq2.keys())
    resvs2.sort()
    n_resvs2 = len(resvs2)

    alg1: List[Optional[str]] = []
    resvs1_ix = 0
    alg1_to_resv: Dict[int, int] = {}
    resv1 = None
    alg2: List[Optional[str]] = []
    resv2 = None
    resvs2_ix = 0
    spatial_ix = 0
    alg2_to_resv: Dict[int, int] = {}

    def next_alg1():
        nonlocal resvs1_ix
        nonlocal alg1
        nonlocal resv1

        if resv1 is not None:
            alg1_to_resv[len(alg1)] = resv1
        alg1.append(next1)
        resvs1_ix += 1

    def skip_alg1():
        nonlocal alg1
        alg1.append(None)

    def next_alg2():
        nonlocal resvs2_ix
        nonlocal alg2

        if resv2 is not None:
            alg2_to_resv[len(alg2)] = resv2
        alg2.append(next2)
        resvs2_ix += 1

    def skip_alg2():
        nonlocal alg2
        alg2.append(None)

    def next_both():
        next_alg1()
        next_alg2()

    def finished():
        return resvs1_ix >= len(resvs1) \
            and resvs2_ix >= len(resvs2)

    while(not finished()):
        resv1 = resvs1[resvs1_ix] if resvs1_ix < n_resvs1 else None
        resv2 = resvs2[resvs2_ix] if resvs2_ix < n_resvs2 else None

        if resv1 is None:
            next1 = None
        else:
            next1 = seq1[resv1]
        
        if resv2 is None:
            next2 = None
        else:
            next2 = seq2[resv2]

        # If one of the sequences has been completed,
        # we just need to append until both sequences
        # are done.
        # It is also the case that if we have already
        # visited all the spatial indexes, then we
        # can just append what remains to the
        # alignment
        if resvs1_ix >= len(resvs1) \
            or resvs2_ix >= len(resvs2) \
            or spatial_ix >= len(spatial_alignment):
            next_both()
            continue

        (ix1, ix2) = spatial_alignment[spatial_ix]
        s_resv1 = int(distances.indexes[0][ix1])
        s_resv2 = int(distances.indexes[1][ix2])

        # If both sequences are at a position matching
        # the spatial alignment, then we append
        # both residues and move forward. We also
        # move to the next spatial index as the current
        # one has been resolved
        if resv1 == s_resv1 and resv2 == s_resv2:
            next_both()
            spatial_ix += 1

        # If both poositions are below the spatial index
        # we need to advance both
        elif (resv1 is not None and resv1 < s_resv1) \
            and (resv2 is not None and resv2 < s_resv2):
            next_both()

        # If one of the sequences has reached a position
        # in the spatial alignment, we advance the other
        # sequence until it also reaches the position in
        # the spatial alignment
        elif resv1 == s_resv1 \
            or (resv2 is not None and resv2 < s_resv2):
            skip_alg1()
            next_alg2()
        elif resv2 == s_resv2 \
            or (resv1 is not None and resv1 < s_resv1):
            next_alg1()
            skip_alg2()
        
        # The last possibility is that the position
        # of both sequences is not part of the spatial
        # alignment. In such case, we advance both
        else:
            next_both()

    with StringIO() as raw_msa:
        raw_msa.write(f">{s1.show()}\n")
        raw_msa.write("".join(v if v is not None else "-" for v in alg1))
        raw_msa.write("\n\n")
        raw_msa.write(f">{s2.show()}\n")
        raw_msa.write("".join(v if v is not None else "-" for v in alg2))
        raw_msa.seek(0)

        msa = AlignIO.read(raw_msa, format='fasta')

    return SpatialAlignment(
        distance_matrix=distances,
        alignment=cast(Any, msa),
        resv_maps=[alg1_to_resv, alg2_to_resv],
        resv_to_seq_maps=[index_by_resv(i) for i in distances.indexes]
    )