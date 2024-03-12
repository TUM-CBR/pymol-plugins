from typing import Any, List, Optional
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, QObject, Qt

from ..data import CavitiesResult, CavityModel

def sort_cavities(cavities: List[CavityModel]):
    new_list = list(cavities)
    new_list.sort(key=lambda cavity: cavity.volume(), reverse=True)

    return new_list

class CavityEntry:
    structure: str
    cavity: CavityModel
    visible: bool

    def __init__(
        self,
        structure: str,
        cavity: CavityModel
    ) -> None:
        self.structure = structure
        self.cavity = cavity
        self.visible = True

class CavityResultsModel(QAbstractTableModel):

    K_STRUCTURE = "structure"
    K_VOLUME = "volume"
    K_DISPLAY = "display"

    COLUMNS = [
        K_DISPLAY,
        K_STRUCTURE,
        K_VOLUME,
    ]

    def __init__(
        self,
        results: CavitiesResult,
        parent: Optional[QObject] = None
    ) -> None:
        super().__init__(parent)
        self.__cavities = [
            CavityEntry(structure, cavity)
            for structure, cavities in results.cavities.items()
            for cavity in sort_cavities(cavities)
        ]

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
    
    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if not index.isValid():
            return
        
        

    

    

