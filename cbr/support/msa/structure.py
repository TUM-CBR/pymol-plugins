from Bio.Align import MultipleSeqAlignment
from typing import List, NamedTuple, Optional

from ...clustal import msa
from ...control import viter
from ...core.pymol import structure
from ...core.pymol.structure import StructureSelection
from ...clustal.Clustal import Clustal, get_clustal
from ...support.msa import from_biopython, Msa

MsaToStructureMap = List[Optional[int]]

def msa_to_structure_position_map(
    sequence_name : str,    
    full_msa : 'Msa | MultipleSeqAlignment',
    structure_sequence : str,
    clustal : Optional[Clustal] = None
)  -> MsaToStructureMap:

    if isinstance(full_msa, MultipleSeqAlignment):
        full_msa = from_biopython(full_msa)

    clustal = clustal or get_clustal()
    sequence = msa.clean_msa_blanks(full_msa[sequence_name])
    result = clustal.run_msa_items(
        [ ("structure", structure_sequence)
        , (sequence_name, sequence)
        ]
    )

    return list(msa.get_relative_positions(full_msa, result))

SequenceToStructureMap = List[int]

def sequence_to_structure_position_map(
    selection : StructureSelection
) -> SequenceToStructureMap:
    return list(sorted(structure.get_pdb_sequence_index(selection).keys()))

class MsaToPymolStructureMap(NamedTuple):
    msa_to_structure : MsaToStructureMap
    sequence_to_structure : SequenceToStructureMap

    def get_pymol_structure_position(self, msa_position: int) -> Optional[int]:

        for structure_pos in viter(self.msa_to_structure[msa_position]):
            return self.sequence_to_structure[structure_pos] + 1

        return None

def msa_to_pymol_structure_map(
    structure_selection: StructureSelection,
    sequence_name : str,    
    full_msa : 'Msa | MultipleSeqAlignment',
    clustal : Optional[Clustal] = None
):
    structure_sequence = structure.get_pdb_sequence(structure_selection)
    msa_to_structure = msa_to_structure_position_map(
        sequence_name,
        full_msa,
        structure_sequence,
        clustal
    )
    sequence_to_structure = sequence_to_structure_position_map(
        structure_selection
    )

    return MsaToPymolStructureMap(
        msa_to_structure = msa_to_structure,
        sequence_to_structure = sequence_to_structure
    )

