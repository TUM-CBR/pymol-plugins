from typing import List, Optional

from ...clustal import msa
from ...core.pymol import structure
from ...core.pymol.structure import StructureSelection
from ...clustal.Clustal import Clustal, get_clustal
from ...support.msa import Msa

def msa_to_structure_position_map(
    sequence_name : str,    
    full_msa : Msa,
    structure_sequence : str,
    clustal : Optional[Clustal] = None
)  -> 'List[int | None]':
    clustal = clustal or get_clustal()
    sequence = msa.clean_msa_blanks(full_msa[sequence_name])
    result = clustal.run_msa_items(
        [ ("structure", structure_sequence)
        , (sequence_name, sequence)
        ]
    )

    return list(msa.get_relative_positions(full_msa, result))

def sequence_to_structure_position_map(
    selection : StructureSelection
) -> List[int]:
    return list(sorted(structure.get_pdb_sequence_index(selection).keys()))