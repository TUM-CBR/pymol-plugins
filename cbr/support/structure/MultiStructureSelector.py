import enum
from PyQt5.QtCore import QByteArray, QMetaType, QModelIndex, QObject, Qt, pyqtSignal, pyqtProperty, pyqtSlot, QAbstractTableModel
from PyQt5.QtWidgets import QComboBox, QDialog, QItemDelegate, QItemEditorCreatorBase, QItemEditorFactory, QWidget
from typing import Any, Dict, List, NamedTuple, Optional, Tuple

from ...core.pymol.structure import StructureSelection
from ...core.pymol import objects

from .Ui_MultiStructureSelector import Ui_MultiStructureSelector

class ValuesSelector(QComboBox):

    def __init__(
        self,
        options: List[str],
        parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(parent)
        self.insertItems(0,options)

    @pyqtProperty(str, user=True)
    def selected_value(self) -> object:
        return self.currentText()

class ValueCreator(QItemEditorCreatorBase):

    def __init__(
        self,
        items: List[str]
    ):
        super().__init__()
        self.__items = items

    def createWidget(self, parent: Optional[QWidget]) -> Optional[QWidget]:
        return ValuesSelector(self.__items, parent=parent)
    
    def valuePropertyName(self) -> QByteArray:
        return QByteArray(b"selected_value")
    
class Factory(QItemEditorFactory):

    def createEditor(self, userType: int, parent: Optional[QWidget]) -> Optional[QWidget]:
        editor = super().createEditor(userType, parent)
        return editor

class StructureAndSequence(NamedTuple):
    structure: StructureSelection
    sequence: str

class StructureEntry:
    structure: StructureSelection
    selected: bool
    sequence: Optional[str]

    def __init__(
        self,
        structure: StructureSelection,
        selected: bool
    ) -> None:
        self.structure = structure
        self.selected = selected
        self.sequence = None

    def set_selected(self, checked: Qt.CheckState) -> None:
        self.selected = checked == Qt.CheckState.Checked

    def set_sequence(self, seq: str) -> None:
        self.sequence = seq

    def checked_state(self) -> Qt.CheckState:
        if self.selected:
            return Qt.CheckState.Checked
        else:
            return Qt.CheckState.Unchecked
        
    def as_structure_and_sequence(self) -> Optional[StructureAndSequence]:
        if self.sequence is None:
            return None
        else:
            return StructureAndSequence(
                self.structure,
                self.sequence
            )

StructureKey = Tuple[str, str, str]

class StructuresTableModel(QAbstractTableModel):

    K_SELECTED = "selected"
    K_NAME = "name"
    K_SEQUENCE = "sequence"
    COLUMNS_BASE = [K_SELECTED, K_NAME]

    selected_values_changed = pyqtSignal()

    @classmethod
    def sequence_column(cls):
        return len(cls.COLUMNS_BASE)

    def __init__(
        self,
        sequences: Optional[List[str]] = None,
        parent: Optional[QObject] = None
    ) -> None:
        super().__init__(parent)
        self.__structures : Dict[StructureKey, StructureEntry] = {}
        self.__keys: List[StructureKey] = list(self.__structures.keys())
        self.__sequences: Optional[List[str]] = sequences
        self.__columns: List[str] = \
            self.COLUMNS_BASE \
                if sequences is None \
                else self.COLUMNS_BASE + [self.K_SEQUENCE]

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
        self.__keys.sort(key = lambda k: self.__structures[k].structure.show())

        self.modelReset.emit()

    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__columns)
    
    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__structures)
    
    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if role != Qt.ItemDataRole.DisplayRole:
            return None
        
        if orientation == Qt.Orientation.Horizontal:
            return self.__columns[section]
        else:
            return None
        
    def flags(self, index: QModelIndex) -> Qt.ItemFlags:

        flags = super().flags(index)
        if index.column() == self.__columns.index(self.K_SELECTED):
            return flags | Qt.ItemFlag.ItemIsUserCheckable
        elif index.column() == self.__columns.index(self.K_SEQUENCE):
            return flags | Qt.ItemFlag.ItemIsEditable
        
        return flags
    
    def __get_entry(self, index: QModelIndex) -> StructureEntry:
        return self.__structures[self.__keys[index.row()]]

    def setData(self, index: QModelIndex, value: Any, role: int = Qt.ItemDataRole.DisplayRole) -> bool:

        column = index.column()
        emit = False
        if role == Qt.ItemDataRole.CheckStateRole and column == self.__columns.index(self.K_SELECTED):
            self.__get_entry(index).set_selected(value)
            emit = True
            
        elif role == Qt.ItemDataRole.EditRole and column == self.__columns.index(self.K_SEQUENCE):
            self.__get_entry(index).set_sequence(value)
            emit = True

        if emit:
            self.dataChanged.emit(index, index)
            self.selected_values_changed.emit()
            return True
        else:
            return super().setData(index, value, role)
    
    def __get_data(self, index: QModelIndex) -> Optional[str]:

        col = index.column()
        value = self.__get_entry(index)

        if col == self.__columns.index(self.K_NAME):
            return value.structure.show()
        if col == self.__columns.index(self.K_SEQUENCE):
            seq = value.sequence
            return seq if seq is not None else "<double click to select>"
        else:
            return None
        
    def __get_checked_state(self, index: QModelIndex) -> Optional[Qt.CheckState]:

        col = index.column()
        value = self.__get_entry(index)

        if col == self.__columns.index(self.K_SELECTED):
            return value.checked_state()
        else:
            return None

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if role == Qt.ItemDataRole.DisplayRole:
            return self.__get_data(index)
        elif role == Qt.ItemDataRole.CheckStateRole:
            return self.__get_checked_state(index)
        else:
            return None
    
    def selected_values(self) -> List[StructureSelection]:
        return [
            entry.structure
            for entry in self.__structures.values() if entry.selected
        ]
    
    def selected_sequences(self) -> List[StructureAndSequence]:
        return [
            value
            for entry in self.__structures.values() if entry.selected
            for value in [entry.as_structure_and_sequence()] if value is not None
        ]

class Mode(enum.Enum):
    Structure = 0
    StructureAndSequence = 1

class MultiStructureSelector(QDialog):

    def __init__(
        self,
        min_structures: Optional[int],
        max_structures: Optional[int],
        parent: Optional[QWidget] = None,
    ) -> None:
        super().__init__(parent)

        self.__ui = Ui_MultiStructureSelector()
        self.__ui.setupUi(self)

        self.__min_structures = min_structures
        self.__max_structures = max_structures
        
        self.__model = StructuresTableModel()
        self.__ui.structuresTable.setModel(self.__model)
        self.__model.selected_values_changed.connect(self.__on_selected_values_changed)

        self.__ui.refreshButton.clicked.connect(self.__on_refresh_clicked)
        self.__ui.cancelButton.clicked.connect(self.__on_cancel)
        self.__ui.acceptButton.clicked.connect(self.__on_accept)
        self.__mode = Mode.Structure

        self.__refresh_structures()
        self.__validate()

    def __set_model(self, model: StructuresTableModel, mode: Mode):
        self.__model.selected_values_changed.disconnect(self.__on_selected_values_changed)

        self.__model = model
        self.__mode = mode
        self.__ui.structuresTable.setModel(self.__model)
        self.__model.selected_values_changed.connect(self.__on_selected_values_changed)

    def exec_with_sequences(self, sequences: List[str]):
        custom_value_delegate = QItemDelegate()
        custom_value_factory = Factory()
        custom_value_factory.registerEditor(QMetaType.Type.UnknownType, ValueCreator(sequences))
        custom_value_delegate.setItemEditorFactory(custom_value_factory)

        self.__ui.structuresTable.setItemDelegateForColumn(StructuresTableModel.sequence_column(), custom_value_delegate)
        self.__set_model(StructuresTableModel(sequences=sequences), Mode.StructureAndSequence)

        self.__refresh_structures()
        return super().exec()

    def exec(self) -> int:
        self.__set_model(StructuresTableModel(), Mode.Structure)
        self.__refresh_structures()
        return super().exec()

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

        if self.__mode == Mode.Structure:
            values = len(self.__model.selected_values())
        else:
            values = len(self.__model.selected_sequences())

        valid = \
            (self.__min_structures is None or values >= self.__min_structures) and \
            (self.__max_structures is None or values <= self.__max_structures)
        
        if valid:
            self.__ui.acceptButton.setEnabled(True)
            self.__ui.infoLabel.setText(None)
        else:
            self.__ui.acceptButton.setEnabled(False)
            self.__ui.infoLabel.setText(self.__validation_text())

    def selected_values(self):
        return self.__model.selected_values()
    
    def selected_sequences(self):
        return self.__model.selected_sequences()