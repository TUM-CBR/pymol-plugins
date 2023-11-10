import math
import pymol

from typing import Dict, List, NamedTuple, Optional, Tuple, TypeVar

class StructureSelection(NamedTuple):
    structure_name : str
    chain_name : Optional[str]
    segment_identifier : Optional[str]

    @property
    def selection(self) -> str:
        selectors = \
            [ self.structure_name
            , self.chain_name and "chain %s" % self.chain_name
            , self.segment_identifier and "segi %s" % self.segment_identifier
            ]

        return " and ".join(s for s in selectors if s)

    def show(self):
        return f"{self.structure_name}/{self.chain_name}/{self.segment_identifier}"

def get_structure_query(structure_name : str, chain : 'str | None' = None) -> str:
    if chain:
        return "(model %s) & (chain %s)" % (structure_name, chain)
    else:
        return "model %s" % structure_name

AnySelection = TypeVar('AnySelection', str, StructureSelection)

def get_selection_sequence_index(selection_obj : AnySelection) -> Dict[int, str]:

    selection = \
        selection_obj \
        if isinstance(selection_obj, str) \
        else selection_obj.selection

    result = []

    pymol.cmd.iterate(
        "%s & guide & alt +A" % selection,
        'result.append((int(resi), oneletter))',
        space={'result': result, 'int': int}
    )
    return dict(result)

def get_pdb_sequence_index(selection : AnySelection) -> Dict[int, str]:
    return get_selection_sequence_index(
        selection
    )

def get_selection_sequece(selection : AnySelection) -> str:
    sequence = get_selection_sequence_index(selection)
    return "".join(
        sequence[k]
        for k in sorted(sequence.keys())
    )

def get_pdb_sequence(selection : StructureSelection) -> str:
    return "".join(
        name
        for buffer in [[]] if pymol.cmd.iterate(
            "(%s) & guide & alt +A & bonded" % selection.selection,
            "buffer.append(oneletter)",
            space={'buffer': buffer}
        ) >= 0
        for name in buffer
    )

AtomPosition = Tuple[float, float, float]

class ResidueDistanceCalculator:

    def __init__(self, selection : str):

        self.__positions : Dict[int, List[AtomPosition]] = dict(
            (k,[])
            for k,_ in get_selection_sequence_index(selection).items()
        )
        def add_atom(resi: int, px: float, py: float, pz: float):
            self.__positions[resi].append((px, py, pz))

        pymol.cmd.iterate_state(
            1,
            selection,
            'add_atom(resv,x,y,z)',
            space={'add_atom': add_atom}
        )

        self.__zero_offset = min(self.__positions.keys())

    def distance(self, resi_a: int, resi_b : int) -> Optional[float]:

        distances = [
            math.sqrt((xa-xb)**2 + (ya-yb)**2 + (za - zb)**2)
            for (xa, ya, za) in self.__positions[resi_a]
            for (xb, yb, zb) in self.__positions[resi_b]
        ]

        # If the resiude is unresolved in the PDB structure, then
        # its position will be unkonwn and the distance cannot be
        # calculated
        if len(distances) > 0:
            return min(distances)
        else:
            return None

    def distance_zero_index(self, resi_a: int, resi_b: int) -> Optional[float]:
        return self.distance(
            resi_a + self.__zero_offset,
            resi_b + self.__zero_offset
        )
