import pymol
from PyQt5.QtCore import Qt, QObject, QModelIndex, pyqtSignal, pyqtSlot, QAbstractTableModel
from PyQt5.QtGui import QColor
from typing import Callable, Dict, Iterable, List, NamedTuple, Optional, Tuple

from ...core import color
from ...core.pymol.structure import StructureSelection

PositionMapping = List[Optional[int]]

class StructureMapping(NamedTuple):
    """
    Type that represents the mapping of a table widget to
    positions in a structure

    Attributes
    ----------
    
    structure : StructureSelection
        The structure for which this mapping is valid
    by_column : Optional[List[int]]
        Represents a column mapping where each column of the table is mapped to
        a position in the structure. The indexes of the list correspond to the
        column numbers and the value to the structure position (resi).
    by_row : Optional[List[int]]
        Represents a mapping where each row of the table is mapped to a
        position in the structure. 
    """
    structure: StructureSelection
    by_column: Optional[PositionMapping] = None
    by_row:  Optional[PositionMapping] = None

class StructureMappings(NamedTuple):
    """
    Data structure representing the mapping of table positions to
    structures.

    Attributes
    ----------

    mappings: List[StructureMappings]
        A list with all of the mappings that have been generated
    """

    mappings: List[StructureMapping]

class StructureEntry(NamedTuple):
    mapping: StructureMapping
    original_colors: Dict[int, int]

BLACK = QColor(0,0,0)

class AbstractStructureTableModel(QAbstractTableModel):

    structure_mappings_signal = pyqtSignal()

    def __init__(self, parent: Optional[QObject]) -> None:
        super().__init__(parent)

        self.structure_mappings_signal.connect(self.__on_structure_mappings_changed)
        self.__current_structures : Optional[List[StructureEntry]] = None

    def __get_structure_mappings__(self) -> Optional[StructureMappings]:
        raise NotImplementedError()

    def __restore_structure(self, structure: StructureEntry):
        structure.mapping.structure.set_colors(structure.original_colors)

    def __restore_current_structures(self):

        current = self.__current_structures
        self.__current_structures = None

        if current is None:
            return
        
        for item in current:
            self.__restore_structure(item)

    def __color_structure(
        self,
        structure: StructureSelection,
        iter: Iterable[Tuple[int, Optional[QColor]]]
    ) -> None:
        for resv, q_color in iter:
            sele = structure.residue_selection([resv])

            if q_color is None:
                q_color = BLACK

            pymol.cmd.color(
                color.to_pymol_color(q_color),
                sele
            )

    def __color_by_column(self, structure: StructureSelection, column_mapping: PositionMapping) -> None:

        self.__color_structure(
            structure,
            (
                (resv, self.__get_structure_color__(None, col))
                for col, resv in enumerate(column_mapping)
                if resv is not None
            )
        )

    def __color_by_row(self, structure: StructureSelection, row_mapping: PositionMapping):

        self.__color_structure(
            structure,
            (
                (resv, self.__get_structure_color__(row, None))
                for row, resv in enumerate(row_mapping)
                if resv is  not None
            )
        )

    def __select_index_by_mapping(
        self,
        selection: StructureSelection,
        indexes: List[QModelIndex],
        mapper: Callable[[QModelIndex], Optional[int]]
    ) -> Optional[StructureSelection]:
        
        if len(indexes) == 0:
            return None


        scope = [
            resv
            for i in indexes
            for resv in [mapper(i)]
            if resv is not None
        ]

        if len(scope) == 0:
            return None

        return selection.scoped(scope)

    def __select_indexes(self, structure: StructureEntry, indexes: List[QModelIndex]) -> Optional[StructureSelection]:

        mapping = structure.mapping

        if mapping.by_column is not None:
            by_column = mapping.by_column
            return self.__select_index_by_mapping(
                mapping.structure,
                indexes,
                lambda i: by_column[i.column()],
            )
        elif mapping.by_row is not None:
            by_row = mapping.by_row
            return self.__select_index_by_mapping(
                mapping.structure,
                indexes,
                lambda i: by_row[i.row()],
            )
        else:
            raise Exception(f"Missing implementations for the mappings definition {mapping}")

    def __apply_structure_coloring(self, structure: StructureEntry):

        coloring = structure.mapping

        if coloring.by_column is not None:
            self.__color_by_column(structure.mapping.structure, coloring.by_column)
        elif coloring.by_row is  not None:
            self.__color_by_row(structure.mapping.structure, coloring.by_row)

    def __iterate_current_structures(self) -> Iterable[StructureEntry]:
        current_structures = self.__current_structures

        if current_structures is None:
            return

        for structure in current_structures:
            yield structure

    def __apply_coloring(self):

        for structure in self.__iterate_current_structures():
            self.__apply_structure_coloring(structure)

    def __to_structure_entry(self, mapping: StructureMapping) -> StructureEntry:
        colors = mapping.structure.get_color_indexes()

        return StructureEntry(
            mapping=mapping,
            original_colors=colors
        )

    def __update_structures(self):
        mappings = self.__get_structure_mappings__()

        if mappings is None:
            self.__current_structures = []
        else:
            self.__current_structures = [
                self.__to_structure_entry(mapping)
                for mapping in mappings.mappings
            ]

    @pyqtSlot()
    def __on_structure_mappings_changed(self):

        self.__restore_current_structures()

        self.__update_structures()
        self.__apply_coloring()
        pass

    def __get_structure_color__(self, row: Optional[int], col: Optional[int]) -> Optional[QColor]:

        if row is None:
            row = 0

        if col is None:
            col = 0

        index = self.index(row, col)
        color = self.data(index, Qt.ItemDataRole.BackgroundColorRole)

        if isinstance(color, QColor):
            return color
        else:
            return None
        
    def update_selection(self, indexes: List[QModelIndex], name: str = "sele"):

        selections = (
            selection.selection
            for entry in self.__iterate_current_structures()   
            for selection in [self.__select_indexes(entry, indexes)]
                if selection is not None
        )

        pymol.cmd.select(
            name,
            " or ".join(selections)
        )