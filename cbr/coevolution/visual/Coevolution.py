from Bio.Align import MultipleSeqAlignment
from typing import Any, List, NamedTuple, Optional, Sequence

from PyQt5.QtCore import QModelIndex, Qt

from ...core.pymol.structure import StructureSelection
from ...core.Context import Context
from ...clustal.Clustal import Clustal
from ...support.structure.AbstractStructureTableModel import AbstractStructureTableModel, StructureMapping, StructureMappings
from ...support.structure.StructuresAlignmentMapper import StructuresAlignmentMapper
from ..data import CoevolutionPosition, CoevolutionResults
from .Ui_coevolution import Ui_Coevolution

class ReferenceStructureEntry(NamedTuple):
    structure: StructureSelection
    reference_sequence_id: str

class CoevolutionResultTable(AbstractStructureTableModel):

    K_RESIDUE_1 = "Local Residue"
    K_RESIDUE_2 = "Other Residue"
    K_SCORE = "Confidence"
    K_SCORE_OCCURENCE = "Occurrence"
    K_SCORE_EXCLUSIVITY = "Exclusivity"
    K_SCORE_CONSERVED = "Conserved"

    COLUMNS = [
        K_RESIDUE_1,
        K_RESIDUE_2,
        K_SCORE,
        K_SCORE_EXCLUSIVITY,
        K_SCORE_CONSERVED
    ]

    def __init__(
        self,
        msa: MultipleSeqAlignment,
        clustal: Clustal
    ):
        self.__msa = msa
        self.__msa_mapper = StructuresAlignmentMapper(msa, clustal)
        self.__positions: List[int] = []
        self.__structures: Optional[Sequence[ReferenceStructureEntry]] = None

    def set_results(self, results: CoevolutionPosition):

        self.__positions = list(results.by_position.keys())
        self.__positions.sort(
            key=lambda pos: results.by_position[pos].score,
            reverse=True
        )

    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__positions)
    
    def columnCount(self, parent: QModelIndex = ...) -> int:
        return len(self.COLUMNS)
    
    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        
        if orientation == Qt.Orientation.Vertical:
            return self.COLUMNS[section]
        else:
            return None
        
    def __create_structure_mapping(self, structure: ReferenceStructureEntry) -> StructureMapping:
        mapping = self.__msa_mapper.map_structure(structure.structure, structure.reference_sequence_id)
        msa_to_resv = {
            pos: resv
            for resv, pos in mapping.resv_to_msa().items()
            if pos is not None
        }

        by_row = [
            msa_to_resv.get(pos)
            for pos in self.__positions
        ]

        return StructureMapping(
            structure.structure,
            by_row = by_row
        )

    def __get_structure_mappings__(self) -> Optional[StructureMappings]:

        if self.__structures is None:
            return None
        
        return StructureMappings(
            mappings = [
                self.__create_structure_mapping(entry)
                for entry in self.__structures
            ]
        )
class Coevolution(Ui_Coevolution):

    def __init__(
        self,
        context: Context
    ) -> None:
        super().__init__()
        self.__ui = Ui_Coevolution()
        self.__ui.setupUi(self)

