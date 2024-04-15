from PyQt5.QtCore import QModelIndex, QObject, Qt, pyqtSignal, pyqtSlot, QAbstractTableModel
from PyQt5.QtWidgets import QDialog, QWidget
from typing import Any, Dict, List, Optional, Tuple

from ...core.pymol.structure import StructureSelection
from ...core.pymol import objects

from .Ui_MultiStructureSelector import Ui_MultiStructureSelector

class StructureEntry:
    structure: StructureSelection
    selected: bool

    def __init__(
        self,
        structure: StructureSelection,
        selected: bool
    ) -> None:
        self.structure = structure
        self.selected = selected

    def set_selected(self, checked: Qt.CheckState) -> None:
        self.selected = checked == Qt.CheckState.Checked

    def checked_state(self) -> Qt.CheckState:
        if self.selected:
            return Qt.CheckState.Checked
        else:
            return Qt.CheckState.Unchecked

StructureKey = Tuple[str, str, str]

class StructuresTableModel(QAbstractTableModel):

    K_SELECTED = "selected"
    K_NAME = "name"
    COLUMNS = [K_SELECTED, K_NAME]

    selected_values_changed = pyqtSignal()

    def __init__(self, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self.__structures : Dict[StructureKey, StructureEntry] = {}
        self.__keys: List[StructureKey] = list(self.__structures.keys())

    def refresh(self):

        not_visited = set(self.__structures.keys())

        for (name, chain, seg) in objects.iter_segis():
            key = (name, chain, seg)

            if key in not_visited:
                not_visited.remove(key)

            if key in self.__structures:
                continue

            self.__structures[key] = StructureEntry(
                structure=StructureSelection(name, chain, seg),
                selected=False
            )

        for item in not_visited:

            if item in self.__structures:
                self.__structures.pop(item)

        self.__keys = list(self.__structures.keys())

        self.modelReset.emit()

    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.COLUMNS)
    
    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__structures)
    
    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if role != Qt.ItemDataRole.DisplayRole:
            return None
        
        if orientation == Qt.Orientation.Horizontal:
            return self.COLUMNS[section]
        else:
            return super().headerData(section, orientation, role)
        
    def flags(self, index: QModelIndex) -> Qt.ItemFlags:

        flags = super().flags(index)
        if index.column() == self.COLUMNS.index(self.K_SELECTED):
            return flags | Qt.ItemFlag.ItemIsUserCheckable
        
        return flags
    
    def __get_entry(self, index: QModelIndex) -> StructureEntry:
        return self.__structures[self.__keys[index.row()]]

    def setData(self, index: QModelIndex, value: Any, role: int = Qt.ItemDataRole.DisplayRole) -> bool:

        if role == Qt.ItemDataRole.CheckStateRole:
            self.__get_entry(index).set_selected(value)
            self.dataChanged.emit(index, index)
            self.selected_values_changed.emit()
            return True

        return super().setData(index, value, role)
    
    def __get_data(self, index: QModelIndex) -> Optional[str]:

        col = index.column()
        value = self.__get_entry(index)

        if col == self.COLUMNS.index(self.K_NAME):
            return value.structure.show()
        else:
            return None
        
    def __get_checked_state(self, index: QModelIndex) -> Optional[Qt.CheckState]:

        col = index.column()
        value = self.__get_entry(index)

        if col == self.COLUMNS.index(self.K_SELECTED):
            return value.checked_state()
        else:
            return None

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if role == Qt.ItemDataRole.DisplayRole:
            return self.__get_data(index)
        elif role == Qt.ItemDataRole.CheckStateRole:
            return self.__get_checked_state(index)

        return super().data(index, role)
    
    def selected_values(self) -> List[StructureSelection]:
        return [
            entry.structure
            for entry in self.__structures.values() if entry.selected
        ]

class MultiStructureSelector(QDialog):

    def __init__(
        self,
        min_structures: Optional[int],
        max_structures: Optional[int],
        parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(parent)

        self.__ui = Ui_MultiStructureSelector()
        self.__ui.setupUi(self)

        self.__min_structures = min_structures
        self.__max_structures = max_structures
        self.__model = StructuresTableModel()
        self.__ui.structuresTable.setModel(self.__model)
        self.__ui.refreshButton.clicked.connect(self.__on_refresh_clicked)
        self.__ui.cancelButton.clicked.connect(self.__on_cancel)
        self.__ui.acceptButton.clicked.connect(self.__on_accept)
        self.__model.selected_values_changed.connect(self.__on_selected_values_changed)

        self.__refresh_structures()
        self.__validate()

    @pyqtSlot()
    def __on_selected_values_changed(self):
        self.__validate()

    @pyqtSlot()
    def __on_refresh_clicked(self):
        self.__refresh_structures()

    @pyqtSlot()
    def __on_accept(self):
        self.accept()

    @pyqtSlot()
    def __on_cancel(self):
        self.reject()

    def __refresh_structures(self):
        self.__model.refresh()

    def __validation_text(self) -> str:

        if self.__min_structures == self.__max_structures:
            return f"You must select {self.__min_structures} structures."
        else:
            low = "0" if self.__min_structures is None else str(self.__min_structures)
            high = "∞" if self.__max_structures is None else str(self.__max_structures)
            return f"You must select between {low} and {high} structures."

    def __validate(self):
        values = self.__model.selected_values()

        valid = \
            (self.__min_structures is None or len(values) >= self.__min_structures) and \
            (self.__max_structures is None or len(values) <= self.__max_structures)
        
        if valid:
            self.__ui.acceptButton.setEnabled(True)
            self.__ui.infoLabel.setText(None)
        else:
            self.__ui.acceptButton.setEnabled(False)
            self.__ui.infoLabel.setText(self.__validation_text())

    def selected_values(self):
        return self.__model.selected_values()