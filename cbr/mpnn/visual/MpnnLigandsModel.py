import pymol
from PyQt5.QtCore import QAbstractTableModel, QItemSelection, QModelIndex, Qt
from typing import Any

from ...core.pymol import objects, structure

class LigandEntry:
    model: str
    chain: str
    selected: bool = False

    def __init__(self, model: str, chain: str, seg: str):
        self.model = model
        self.chain = chain
        self.seg = seg
        self.selected = False

    @classmethod
    def get_current_models(cls):

        for item in objects.iter_chains():
            (name, chain) = item
            for seg in item.iter_segs():
                yield LigandEntry(name, chain, seg)

    def as_selection(self) -> structure.StructureSelection:
        return structure.StructureSelection(self.model, self.chain, self.seg)

    def show(self):
        return "%s/%s/%s" % (self.model, self.chain, self.seg)

class MpnnLigandModel(QAbstractTableModel):
    COLUMN_SELECTED = "Selected"
    COLUMN_OBJECT = "Object"

    COLUMN_NAMES = [COLUMN_SELECTED, COLUMN_OBJECT]

    def __init__(self):
        super(MpnnLigandModel, self).__init__()
        self.__ligands = list(LigandEntry.get_current_models())


    def rowCount(self, parent):
        return len(self.__ligands)

    def columnCount(self, parent):
        return len(self.COLUMN_NAMES)

    def flags(self, index):
        flags = super(MpnnLigandModel, self).flags(index)
        if index.column() == 0:
            return Qt.ItemIsUserCheckable | flags
        else:
            return flags

    def setData(self, index: QModelIndex, value: Any, role: Qt.ItemDataRole = Qt.DisplayRole) -> bool:
        if role == Qt.CheckStateRole and index.column() == 0:
            self.__ligands[index.row()].selected = value
            self.dataChanged.emit(index, index)
            return True
        else:
            return super(MpnnLigandModel, self).setItemData(index, value)

    def data_checked_state_role(self, index):
        if self.COLUMN_NAMES[index.column()] == self.COLUMN_SELECTED:
            ligand = self.__ligands[index.row()]
            return Qt.Checked if ligand.selected else Qt.Unchecked
        else:
            return None

    def data_display_role(self, index):
        if self.COLUMN_NAMES[index.column()] == self.COLUMN_OBJECT:
            ligand = self.__ligands[index.row()]
            return ligand.show()
        else:
            return None

    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            return self.data_display_role(index)
        elif role == Qt.CheckStateRole:
            return self.data_checked_state_role(index)

    def refresh(self):
        self.beginResetModel()
        self.__ligands = list(LigandEntry.get_current_models())
        self.endResetModel()

    def get_selected_ligands(self):
        return [ligand.as_selection() for ligand in self.__ligands if ligand.selected]

    def on_selection(self, selected: QItemSelection):

        selection_list = [
            f"({item.as_selection().selection_any})"
            for index in selected.indexes()
            for item in [self.__ligands[index.row()]]
                if index.column() == 0
        ]

        selection = " or ".join(selection_list)
        
        pymol.cmd.select("mpnn_ligand", selection)

