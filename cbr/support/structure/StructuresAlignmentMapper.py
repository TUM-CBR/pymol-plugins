from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from typing import Dict, NamedTuple, Optional, Tuple

from ...clustal.Clustal import Clustal
from ...core.pymol.structure import StructureSelection

K_STRUCTURE_ID = "structure"

class StructureResidue(NamedTuple):
    resv: int
    resv_name: str

class StructureAlignmentEntry(NamedTuple):
    structure: StructureSelection
    reference_sequence_id: str
    structure_sequence: Dict[int, Tuple[int, str]]
    mapper_msa: MultipleSeqAlignment
    target_msa: MultipleSeqAlignment

    def resv_to_msa(self) -> Dict[int, Optional[int]]:

        mapper_structure_id = self.mapper_msa[1].id
        structure_to_mapper = self.mapper_msa.alignment.inverse_indices[0]
        mapper_to_reference_seq = self.mapper_msa.alignment.indices[1]

        mapper_seq_ix = next(
            (ix for ix,seq in enumerate(self.target_msa) if seq.id == mapper_structure_id)
        )
        mapper_seq_mappings = self.target_msa.alignment.inverse_indices[mapper_seq_ix]

        return {
            resv: None if ref_seq_pos < 0 else int(mapper_seq_mappings[ref_seq_pos])
            for resv, (pos, _) in self.structure_sequence.items()
            for ref_seq_pos in [mapper_to_reference_seq[structure_to_mapper[pos]]]
        }
    
    def msa_to_resv(self) -> Dict[int, StructureResidue]:
        return {
            msa: StructureResidue(resv, self.structure_sequence[resv][1])
            for resv, msa in self.resv_to_msa().items()
            if msa is not None
        }

class StructuresAlignmentMapper(NamedTuple):

    msa: MultipleSeqAlignment
    clustal: Clustal

    def map_structure(
        self,
        structure: StructureSelection,
        reference_sequence_id: str
    ):
        
        reference_sequence: Optional[SeqRecord] = next(
            (seq for seq in self.msa if seq.id == reference_sequence_id),
            None
        )

        if reference_sequence is None:
            raise ValueError(f"The sequence {reference_sequence_id} is not in the alignment.")
        
        seq_dict = structure.get_sequence()
        seq_resv = list(seq_dict.keys())
        seq_resv.sort()

        structure_seq = "".join(seq_dict[i] for i in seq_resv)
        structure_record = SeqRecord(Seq(structure_seq), id=K_STRUCTURE_ID)

        reference_record = SeqRecord(
            reference_sequence._seq.replace("-", ""),
            reference_sequence.id
        )

        mapper_msa = self.clustal.run_msa_seqs([structure_record, reference_record])

        entry = StructureAlignmentEntry(
            structure=structure,
            reference_sequence_id=reference_sequence_id,
            structure_sequence = {
                resv: (pos, seq_dict[resv])
                for pos, resv in enumerate(seq_resv)
            },
            mapper_msa=mapper_msa,
            target_msa=self.msa
        )

        return entry



