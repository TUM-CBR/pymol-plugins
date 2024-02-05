from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, QUrl, Qt, pyqtBoundSignal, pyqtSignal, pyqtSlot
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtWidgets import QTableWidgetItem, QWidget
from pymol import cmd
from tempfile import TemporaryDirectory
from typing import Any, cast, List, Optional, Protocol

from ...core import assertions
from ...core.Context import Context
from ...core.visual import StructureSelection
from ..support import parsers
from ..support.parsers import ChimeraId, SchemaEnergy
from .Ui_SchemaEnergyViewer import Ui_SchemaEnergyViewer

class ChimeraMapper(Protocol):

    def mappings_changed(self) -> Optional[pyqtBoundSignal]:
        raise NotImplementedError()
    
    def map_chimera(self, chimera: ChimeraId) -> str:
        raise NotImplementedError()

class ChimerasModel(QAbstractTableModel):

    class DefaultChimeraMapper:

        def mappings_changed(self):
            return None
        
        def map_chimera(self, chimera: ChimeraId) -> str:
            return "".join(str(i) for i in chimera)

    CHIMERA_COLUMN = 0
    ENERGY_COLUMN = 1
    MUTATIONS_COLUMN = 2
    COLUMN_COUNT = MUTATIONS_COLUMN + 1

    def __init__(
        self,
        results : parsers.SchemaEnergy,
        chimera_name_mapping: Optional[ChimeraMapper] = None
    ) -> None:
        super().__init__()

        self.__results = results
        self.__chimera_name_mapping = chimera_name_mapping or self.DefaultChimeraMapper()
        on_mappings_changed = self.__chimera_name_mapping.mappings_changed()

        if on_mappings_changed is not None:
            on_mappings_changed.connect(self.modelReset)

        average_entry = SchemaEnergy.Entry(
            chimera = None,
            energy = results.average_energy,
            mutations = results.average_mutations,
            extra = {}
        )
        self.__entries = [average_entry] + list(self.__results.entries)
        self.__extra_columns = list(
            set(
                col
                for entry in self.__entries
                for col in entry.extra 
            )
        )
        self.__all_columns = ["Chimera", "Energy", "Mutations"] + list(self.__extra_columns)

    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            return self.__all_columns[section]
        
        return super().headerData(section, orientation, role)
    
    def columnCount(self, parent: Any = None) -> int:
        return len(self.__all_columns)
    
    def rowCount(self, parent: Any = None) -> int:
        return len(self.__entries)
    
    def __get_entry(self, row: int):
        return self.__entries[row]

    def __display_chimera(self, chimera: Optional[List[int]]) -> str:
        if chimera is None:
            return "<Average>"
        
        return self.__chimera_name_mapping.map_chimera(chimera)

    def __get_display_data_by_ix(
        self,
        row: int,
        column: int
    ) -> Optional[Any]:
        entry = self.__get_entry(row)

        if column == self.CHIMERA_COLUMN:
            return self.__display_chimera(entry.chimera)
        elif column == self.ENERGY_COLUMN:
            return entry.energy
        elif column == self.MUTATIONS_COLUMN:
            return entry.mutations
        else:
            field = self.__all_columns[column]
            return entry.extra.get(field)
    
    def __get_display_data(self, index: QModelIndex):
        return self.__get_display_data_by_ix(
            index.row(),
            index.column()
        )

    def data(
        self,
        index: QModelIndex,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if not index.isValid():
            return None
        
        if role == Qt.ItemDataRole.DisplayRole:
            return self.__get_display_data(index)
        
        return None
    
class ParentsNamesModel(QAbstractTableModel):

    mappings_changed_signal = pyqtSignal()

    def __init__(
        self,
        parents: MultipleSeqAlignment,
        results: SchemaEnergy
    ) -> None:
        super().__init__()

        self.__parents = parents
        self.__chimeras_count = assertions.assert_same_length(*
            (entry.chimera for entry in results.entries if entry.chimera is not None)
        )
        self.__mappings = [
            [str(i+1) for _1 in range(0, self.__chimeras_count)]
            for i,_2 in enumerate(self.__parents)
        ]
        self.__headers = [seq.id for seq in parents]

    def rowCount(self, parent: Any = None) -> int:
        return len(self.__parents)
    
    def columnCount(self, parent: Any = None) -> int:
        return self.__chimeras_count
    
    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if orientation == Qt.Orientation.Vertical and role == Qt.ItemDataRole.DisplayRole:
            return self.__headers[section]

        return super().headerData(section, orientation, role)

    def flags(self, index: QModelIndex) -> Qt.ItemFlags:
        return super().flags(index) | Qt.ItemFlag.ItemIsEditable
    
    def setData(
        self,
        index: QModelIndex,
        value: Any,
        role: int = Qt.ItemDataRole.EditRole
    ) -> bool:
        row = index.row()
        col = index.column()
        self.__mappings[row][col] = value

        self.dataChanged.emit(index, index)
        self.mappings_changed_signal.emit()
        return True
    
    def map_chimera(self, chimera: List[int]) -> str:

        mappings = [self.__mappings[chimera - 1][pos] for pos, chimera in enumerate(chimera)]
        return "".join(mappings)
    
    def mappings_changed(self):
        return self.mappings_changed_signal
    
    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if not index.isValid():
            return None
        
        if role == Qt.ItemDataRole.DisplayRole:
            row = index.row()
            col = index.column()

            return self.__mappings[row][col]
        
        return None

class SchemaEnergyViewer(QWidget):

    CONTACTS_HEADERS = [
        "i",
        "j",
        "pdb i",
        "pdb j",
    ]

    def __init__(
        self,
        context : Context,
        structure_seleciton : StructureSelection,
        raw_results : str,
        raw_contacts : str,
        raw_parents_msa: str,
        results_folder : TemporaryDirectory[Any],
        *args : Any,
        **kwargs : Any
    ):
        super(SchemaEnergyViewer, self).__init__(*args, **kwargs)
        self.__ui = Ui_SchemaEnergyViewer()
        self.__ui.setupUi(self)
        self.__results = parsers.parse_schema_energy(raw_results)
        self.__names_model = ParentsNamesModel(
            cast(MultipleSeqAlignment, AlignIO.read(raw_parents_msa, format="clustal")),
            self.__results
        )

        result_model = ChimerasModel(
            self.__results,
            self.__names_model
        )
        self.__ui.resultsTable.setModel(result_model)
        self.__ui.chimeraNamesTable.setModel(self.__names_model)
        self.__contacts = parsers.parse_schema_contacts(raw_contacts)
        self.__structure_selection = structure_seleciton
        self.__render_everything()
        self.__results_folder = results_folder

    def __del__(self):
        self.__results_folder.cleanup()

    @pyqtSlot()
    def on_openResultsFolder_clicked(self):
        QDesktopServices.openUrl(QUrl.fromLocalFile(self.__results_folder.name))

    @pyqtSlot()
    def on_contactsTable_itemSelectionChanged(self):

        selected = self.__ui.contactsTable.selectedIndexes()
        contacts = self.__contacts
        residues = " | ".join(
            "resi %i | resi %i" % (contacts.contacts[sel.row()].pdb_i, contacts.contacts[sel.row()].pdb_j)
            for sel in selected
        )

        if(len(residues) > 0):
            cmd.select(
                'SCHEMA_CONTACTS',
                "%s & (%s)" % (self.__structure_selection.selection, residues)
            )

    def __render_everything(self):

        interactions_meta = self.__contacts.interactions
        contacts_table = self.__ui.contactsTable
        contacts_table.reset()
        contacts = self.__contacts
        contacts_table.setRowCount(len(contacts.contacts))
        contacts_table.setColumnCount(4 + len(interactions_meta))
        contacts_table.setHorizontalHeaderLabels(
            SchemaEnergyViewer.CONTACTS_HEADERS + \
            [interaction.name for interaction in interactions_meta]
        )
        residues = {}
        cmd.iterate(
            self.__structure_selection.selection,
            'residues[to_int(resi)] = resn',
            space = { 'residues' : residues, 'to_int': int }
        )

        for (i, contact) in enumerate(contacts.contacts):
            contacts_table.setItem(i, 0, QTableWidgetItem(str(contact.seq_i)))
            contacts_table.setItem(i, 1, QTableWidgetItem(str(contact.seq_j)))
            
            pdb_i = contact.pdb_i
            pdb_i_name = residues.get(pdb_i) or "?"
            contacts_table.setItem(
                i,
                2,
                QTableWidgetItem("%i (%s)" % (pdb_i, pdb_i_name))
            )

            pdb_j = contact.pdb_j
            pdb_j_name = residues.get(pdb_j) or "?"
            contacts_table.setItem(
                i,
                3,
                QTableWidgetItem("%i (%s)" % (pdb_j, pdb_j_name))
            )

            interactions_summary = dict((meta.name, 0.0) for meta in interactions_meta)
            for interaction in contact.interactions:
                interactions_summary[interaction.name] += interaction.strength

            for (j, meta) in enumerate(interactions_meta):
                value = round(interactions_summary[meta.name], ndigits=2)
                contacts_table.setItem(i, j + 4, QTableWidgetItem(str(value)))


