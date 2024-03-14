from pymol import cmd
from typing import Any, Dict, List, Optional, Sequence
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, QObject, Qt

from ...core.pymol.structure import StructureSelection
from ..data import CavitiesResult, CavityModel

def sort_cavities(cavities: List[CavityModel]):
    new_list = list(enumerate(cavities))
    new_list.sort(key=lambda item: item[1].volume(), reverse=True)

    return new_list

class CavityEntry:
    name: str
    index: int
    structure: StructureSelection
    cavity: CavityModel
    visible: bool

    def __init__(
        self,
        name: str,
        index: int,
        structure: StructureSelection,
        cavity: CavityModel
    ) -> None:
        self.name = name
        self.index = index
        self.structure = structure
        self.cavity = cavity
        self.visible = True

    @property
    def selection(self) -> str:
        return f"model {self.name} and chain A{self.index}"

    def surrounding_selection(self, distance: float) -> str:
        distance = round(distance, 1)
        selection = f"({self.structure.selection}) within {distance} of ({self.selection})"
        return selection
    
    def set_visible(self, visible: bool):
        self.visible = visible

        if visible == True:
            cmd.show(
                representation = 'mesh',
                selection = self.selection
            )
        else:
            cmd.hide(
                selection = self.selection
            )

class CavityResultsModel(QAbstractTableModel):

    K_STRUCTURE = "structure"
    K_VOLUME = "volume (Å³)"
    K_DISPLAY = "display"

    COLUMNS = [
        K_DISPLAY,
        K_STRUCTURE,
        K_VOLUME,
    ]

    def __init__(
        self,
        results: CavitiesResult,
        structures: Dict[str, StructureSelection],
        parent: Optional[QObject] = None
    ) -> None:
        super().__init__(parent)
        self.__cavities = [
            CavityEntry(
                results.get_cavity_name(structure),
                i,
                structures[structure],
                cavity,
            )
            for structure, cavities in results.cavities.items()
            for i,cavity in sort_cavities(cavities)
        ]

    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if role == Qt.ItemDataRole.DisplayRole:
            if orientation == Qt.Orientation.Horizontal:
                return self.COLUMNS[section]

        return super().headerData(section, orientation, role)

    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__cavities)
    
    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.COLUMNS)
    
    def __column(self, index: QModelIndex):
        return self.COLUMNS[index.column()]

    def flags(self, index: QModelIndex) -> Qt.ItemFlags:

        flags = super().flags(index)

        if self.__column(index) == self.K_DISPLAY:
            return flags | Qt.ItemFlag.ItemIsUserCheckable

        return flags
    
    def __get_record(self, index: QModelIndex) -> CavityEntry:
        return self.__cavities[index.row()]

    def __get_column(self, index: QModelIndex) -> str:
        return self.COLUMNS[index.column()]

    def __display_role(self, index: QModelIndex) -> Any:

        record = self.__get_record(index)
        column = self.__get_column(index)

        if column == self.K_STRUCTURE:
            return record.structure.show()
        elif column == self.K_VOLUME:
            return record.cavity.volume()

        return None
    
    def __check_state_role(self, index: QModelIndex):
        record = self.__get_record(index)
        column = self.__get_column(index)

        if column == self.K_DISPLAY:
            return Qt.CheckState.Checked if record.visible else Qt.CheckState.Unchecked
        
        return None
    
    def update_visibility(self, visibility: bool) -> None:

        for record in self.__cavities:
            record.set_visible(visibility)

        self.dataChanged.emit(
            self.index(0,0),
            self.index(self.rowCount(), self.columnCount())
        )
    
    def setData(self, index: QModelIndex, value: Any, role: int = Qt.ItemDataRole.EditRole) -> bool:

        if role == Qt.ItemDataRole.CheckStateRole:
            record = self.__get_record(index)
            record.set_visible(True if value == Qt.CheckState.Checked else False)
            self.dataChanged.emit(index, index)
            return True

        return super().setData(index, value, role)
    
    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if not index.isValid():
            return
        
        if role == Qt.ItemDataRole.DisplayRole:
            return self.__display_role(index)
        if role == Qt.ItemDataRole.CheckStateRole:
            return self.__check_state_role(index)
        
    def items_selected(self, indexes: Sequence[QModelIndex], distance: float):

        selections = " or ".join([
            f"({self.__get_record(index).surrounding_selection(distance)})"
            for index in indexes
        ])

        cmd.select(
            "sele",
            f"byres ({selections})"
        )