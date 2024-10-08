import math
import pymol

from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple, TypeVar

class ColorIndex(NamedTuple):
    resv: int
    name: str
    chain: str

StructureColors = Dict[ColorIndex, int]

class StructureSelection(NamedTuple):
    structure_name : str
    chain_name : Optional[str]
    segment_identifier : Optional[str]
    residues : Optional[List[int]] = None

    def __residues_query(self, resis: Iterable[int]) -> str:
        return " or ".join(
            f"resi {i}" if i >= 0 else f"resi \\-{abs(i)}"
            for i in resis
        )

    @property
    def base_query(self) -> str:
        return self.base_query_any()

    def base_query_any(self, modifiers: str = " and polymer") -> str:
        selectors = \
            [ f"model {self.structure_name}"
            , self.chain_name and "chain %s" % self.chain_name
            , self.segment_identifier and "segi %s" % self.segment_identifier
            ]
        
        if self.residues is not None:
            selectors += [f"({self.__residues_query(self.residues)})"]

        items = " and ".join(s for s in selectors if s)
        return f"{items} {modifiers}"
    
    @property
    def selection(self) -> str:
        return f"byres ({self.base_query})"

    @property
    def selection_any(self) -> str:
        return self.base_query_any(modifiers="")
    
    def scoped(self, residues: Sequence[int]) -> 'StructureSelection':
        existing = self. residues

        if existing is not None:
            scope = [r for r in existing if r in residues]
        else:
            scope = residues

        return self._replace(residues = scope)

    def show(self):
        return "/".join(
            part
            for part in [self.structure_name, self.chain_name, self.segment_identifier]
                if part is not None
        )
    
    def residue_selection(self, resis: Iterable[int]) -> str:
        resi_selection = self.__residues_query(resis)

        return f"byres ({self.base_query} and ({resi_selection}))"
    
    def get_sequence(self) -> Dict[int, str]:

        query = f"{self.base_query} and name CA"
        result: Dict[int, str] = {}

        def apply(resv: int, res: str):
            result[resv] = res

        pymol.cmd.iterate(
            query,
            'apply(resv, oneletter)',
            space={'apply': apply}
        )

        return result
    
    def get_color_indexes(self) -> Dict[ColorIndex, int]:
        colors : Dict[ColorIndex, int] = {}

        def apply(resv: int, name: str, chain: str, color: int):
            colors[ColorIndex(resv, name, chain)] = color

        pymol.cmd.iterate(
            self.base_query,
            'apply(resv, name, chain, color)',
            space={'apply': apply}
        )

        return colors
    
    def set_colors(self, color_indexes: Dict[ColorIndex, int]):

        def apply(resv: int, name: str, chain: str, color: int):
            return color_indexes.get(ColorIndex(resv, name, chain), color)

        pymol.cmd.alter(
            self.base_query,
            'color = apply(resv, name, chain, color)',
            space = {'apply': apply}
        )

        pymol.cmd.recolor(self.base_query)

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

    result: List[Tuple[int, str]] = []

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

def get_pdb_dominant_color_index(
    selection : StructureSelection
) -> int:
    colors : Dict[int, int] = dict()

    def consume(color: int):
        if color not in colors:
            colors[color] = 0

        colors[color] += 1

    pymol.cmd.iterate(
        selection.selection,
        'consume(color)',
        space={'consume': consume}
    )

    favorite_ix = max(colors.items(), key=lambda kv: kv[1])[0]
    return favorite_ix

def get_pdb_sequence(
    selection : StructureSelection,
    include_non_bonded: bool = False
) -> str:

    parts = [
        f"({selection.selection})",
        "guide",
        "alt +A"
    ]

    if not include_non_bonded:
        parts.append("bonded")

    return "".join(
        name
        for buffer in [[]] if pymol.cmd.iterate(
            " & ".join(parts),
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
