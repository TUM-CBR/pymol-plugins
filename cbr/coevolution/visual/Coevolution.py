from Bio.Align import MultipleSeqAlignment
import numpy as np
from typing import Any, Dict, List, NamedTuple, Optional, Sequence

from PyQt5.QtCore import QModelIndex, Qt
from PyQt5.QtGui import QColor

from ...core.pymol.structure import StructureSelection
from ...core.Context import Context
from ...clustal.Clustal import Clustal
from ...support.structure.AbstractStructureTableModel import AbstractStructureTableModel, StructureMapping, StructureMappings
from ...support.display.sequence import RESIDUE_COLORS
from ...support.structure.StructuresAlignmentMapper import StructuresAlignmentMapper
from ..data import CoevolutionEntry, CoevolutionPosition, CoevolutionResults
from .Ui_coevolution import Ui_Coevolution

class ReferenceStructureEntry(NamedTuple):
    structure: StructureSelection
    reference_sequence_id: str
    msa_to_resv: Dict[int, int]

    @classmethod
    def new(
        cls,
        structure: StructureSelection,
        reference_sequence_id: str,
        msa_mapper: StructuresAlignmentMapper
    ) -> 'ReferenceStructureEntry':
        mapping = msa_mapper.map_structure(structure, reference_sequence_id)
        msa_to_resv = {
            pos: resv
            for resv, pos in mapping.resv_to_msa().items()
            if pos is not None
        }

        return ReferenceStructureEntry(
            structure,
            reference_sequence_id,
            msa_to_resv
        )

class CoevolutionResultTable(AbstractStructureTableModel):

    K_POSITION = "MSA Position"
    K_RESIDUE_1 = "Local Residue"
    K_RESIDUE_2 = "Other Residue"
    K_SCORE = "Confidence"
    K_SCORE_OCCURENCE = "Occurrence"
    K_SCORE_EXCLUSIVITY = "Exclusivity"
    K_SCORE_CONSERVED = "Conserved"

    SCORE_COLOR_LOW = np.array([128,128,0], dtype=np.int8)
    SCORE_COLOR_HIGH = np.array([0, 255, 0], dtype=np.int8)

    COLUMNS = [
        K_POSITION,
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
        self.__msa_mapper = StructuresAlignmentMapper(msa, clustal)
        self.__positions: List[int] = []
        self.__structures: Optional[Sequence[ReferenceStructureEntry]] = None
        self.__results: Optional[CoevolutionPosition] = None
        self.__structure_positions_columns: List[str] = []

    def columns(self) -> Sequence[str]:
        return self.__structure_positions_columns + self.COLUMNS

    def set_results(self, results: CoevolutionPosition):

        self.__positions = list(results.by_position.keys())
        self.__positions.sort(
            key=lambda pos: results.by_position[pos].score,
            reverse=True
        )
        self.__results = results

    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__positions)
    
    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.columns())
    
    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        
        if orientation == Qt.Orientation.Vertical:
            return self.columns()[section]
        else:
            return None
        
    def __create_structure_mapping(self, structure: ReferenceStructureEntry) -> StructureMapping:
        msa_to_resv = structure.msa_to_resv

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
        
        mappings = StructureMappings(
            mappings = [
                self.__create_structure_mapping(entry)
                for entry in self.__structures
            ]
        )

        self.__structure_positions_columns = [
            f"Position ({mapping.structure.show()})"
            for mapping in mappings.mappings
        ]

        return mappings
    
    def __get_position(self, index: QModelIndex) -> Optional[int]:
        row = index.row()

        if len(self.__positions) <= row:
            return None
        
        return self.__positions[row]
    
    def __get_record(self, index: QModelIndex) -> Optional[CoevolutionEntry]:

        column = index.column()

        if column < len(self.__structure_positions_columns):
            return None

        results = self.__results
        key = self.__get_position(index)

        if key is None or results is None:
            return None

        return results.by_position[key]
    
    def __get_structure_entry(self, index: QModelIndex) -> Optional[ReferenceStructureEntry]:

        column = index.column()
        structures = self.__structures

        if column >= len(self.__structure_positions_columns) or structures is None:
            return None
        
        return structures[column]
    
    def __get_structure_position(self, index: QModelIndex) -> Optional[int]:

        structure_entry = self.__get_structure_entry(index)

        if structure_entry is None:
            return None

        position = self.__get_position(index)

        return None if position is None else structure_entry.msa_to_resv[position]
    
    def __data_display_role(self, index: QModelIndex) -> Any:

        structure_position = self.__get_structure_position(index)

        if structure_position is not None:
            return structure_position

        column = self.columns()[index.column()]

        if column == self.K_POSITION:
            return self.__get_position(index)

        record = self.__get_record(index)

        if record is None:
            return None

        if column == self.K_RESIDUE_1:
            return record.residue_1
        elif column == self.K_RESIDUE_2:
            return record.residue_2
        elif column == self.K_SCORE:
            return record.score
        elif column == self.K_SCORE_OCCURENCE:
            return record.score_occurence
        elif column == self.K_SCORE_EXCLUSIVITY:
            return record.score_exclusivity
        elif column == self.K_SCORE_CONSERVED:
            return record.score_conserved
        else:
            return None
    
    def __get_score_color(self, record: CoevolutionEntry) -> Optional[QColor]:

        color_vector = (
            self.SCORE_COLOR_LOW + \
            ((self.SCORE_COLOR_HIGH - self.SCORE_COLOR_LOW) * record.score)
        ).astype(np.int8)

        return QColor(
            r = color_vector[0],
            g = color_vector[1],
            b = color_vector[2]
        )
    
    def __get_structure_color__(self, row: Optional[int], col: Optional[int]) -> Optional[QColor]:

        if row is None:
            return None
        
        index = self.index(row, self.columns().index(self.K_SCORE))
        record = self.__get_record(index)

        if record is None:
            return None

        return self.__get_score_color(record)
    
    def __get_residue_color(self, residue: str)  -> Optional[QColor]:

        return RESIDUE_COLORS.get(residue.upper())

    def __data_background_color_role(self, index: QModelIndex) -> Any:

        column = self.columns()[index.column()]
        record = self.__get_record(index)

        if record is None:
            return None

        if column == self.K_SCORE:
            return self.__get_score_color(index)
        elif column == self.K_RESIDUE_1:
            return self.__get_residue_color(record.residue_1)
        elif column == self.K_RESIDUE_2:
            return self.__get_residue_color(record.residue_2)
        else:
            return None

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if role == Qt.ItemDataRole.DisplayRole:
            return self.__data_display_role(index)
        elif role == Qt.ItemDataRole.BackgroundColorRole:
            return self.__data_background_color_role(index)
        else:
            return None

class Coevolution(Ui_Coevolution):

    def __init__(
        self,
        context: Context
    ) -> None:
        super().__init__()
        self.__ui = Ui_Coevolution()
        self.__ui.setupUi(self)

