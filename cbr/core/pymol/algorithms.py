import numpy as np
from numpy import linalg
from numpy.typing import NDArray
from pymol import cmd
from typing import Any, Generic, List, TypeVar

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

def cx_centroids(
    s: StructureSelection,
    state: int = 1
) -> StructureVector[np.float64]:
    
    resv_maping:List[int] = []
    ca_coords: List[List[List[float]]] = []

    def apply(resv: int, x: float, y: float, z: float):

        if len(resv_maping) == 0 or resv_maping[-1] != resv:

            while len(ca_coords) > 0 and len(ca_coords[-1]) < 2:
                ca_coords[-1].append([0,0,0])

            resv_maping.append(resv)
            ca_coords.append([])

        if len(ca_coords[-1]) < 2:
            ca_coords[-1].append([x, y, z])

    cmd.iterate_state(
        state,
        f"{s.base_query} & (name CA or name CB)",
        'apply(resv, x, y, z)',
        space={'apply': apply}
    )

    return StructureVector(
        selections=[s],
        indexes=[np.array(resv_maping)],
        array=np.array(ca_coords)
    )

def alpha_distance_matrix(
    s1: StructureSelection,
    s2: StructureSelection,
    state1: int = 1,
    state2: int = 1
) -> StructureVector[np.float64]:
    c1 = cx_centroids(s1, state1)
    c2 = cx_centroids(s2, state2)
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