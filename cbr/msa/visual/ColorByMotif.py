import pymol
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, pyqtSlot, Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QFileDialog, QWidget
import re
from typing import Dict, Iterable, Iterator, List, NamedTuple, Optional

from ...control import viter
from ...core import color
from ...core.color import RgbColor
from ...core.Qt.QtWidgets import export_table_view_to_csv, with_error_handler
from ...support.msa import Msa
from ...core.pymol.visual.PymolResidueSelectionModel import pymol_residue_selection_model

from .MsaContext import MsaContext
from .support import msa_to_pymol_structure_map, MsaToPymolStructureMap
from .Ui_ColorByMotif import Ui_ColorByMotif

RESIDUE_ITEM_DATA_ROLE = Qt.UserRole

class MotifMatch(NamedTuple):
    start : int
    end : int

class MotifMatches(NamedTuple):
    matches : List[MotifMatch]

    def add_match(self, match : re.Match):
        self.matches.append(
            MotifMatch(
                start=match.start(),
                end = match.end()
            )
        )

class SequenceMatchEntry(NamedTuple):
    sequence_name : str
    sequence : str
    matches : Dict[str, MotifMatches]

    def add_match(self, motif : str, match : re.Match):

        if motif not in self.matches:
            self.matches[motif] = MotifMatches(list())

        self.matches[motif].add_match(match)

    def count(self, motif : str) -> int:

        match = self.matches.get(motif)

        if match is None:
            return 0

        return len(match.matches)

def get_color_rgb(index: int) -> color.RgbColor:
    return color.get_color_rgb(index, options=color.distinct_colors_5)

def get_qt_color(index: int) -> QColor:
    return QColor(*color.get_qt_color(index, options=color.distinct_colors_5))

class PatternEntry(NamedTuple):
    pattern : re.Pattern
    color_index : int

    @property
    def rgb_color(self) -> RgbColor:
        return get_color_rgb(self.color_index)

    @property
    def qt_color(self) -> QColor:
        return get_qt_color(self.color_index)

    def finditer(self, input : str) -> Iterator[re.Match]:
        return self.pattern.finditer(input)

MotifsDict = Dict[str, PatternEntry]

class MatchMotifsResult(NamedTuple):
    by_msa_position_results : Dict[int, Dict[str, int]]
    by_sequence_name_results : Dict[str, SequenceMatchEntry]
    motifs : MotifsDict

MOTIF_ID_ROLE = Qt.UserRole

KNOWN_MOTIFS : Dict[str, re.Pattern] = {
    'glycosolation': re.compile(r"n\-*[a-z]\-*(t|s)", re.IGNORECASE)
}

VALID_MOTIF_NAME = re.compile(r"\w(\w|\s)*")

class MotifsModel(QAbstractTableModel):

    def __init__(
        self,
        selection : MotifsDict
    ):
        super().__init__()

        # Copy the dictionary, don't trust the outside
        # forces not to modify it
        self.__selection = dict(selection)
        self.__keys = list(self.__selection.keys())

    def rowCount(self, parent = None) -> int:
        return len(self.__keys)

    def columnCount(self, parent = None) -> int:
        return 2

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            headers = ["Name", "Pattern"]
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def data(self, index: QModelIndex, role=Qt.DisplayRole):

        if not index.isValid():
            return None

        key = self.__keys[index.row()]
        if role == Qt.DisplayRole:
            cols : List[str] = [
                key,
                self.__selection[key].pattern.pattern
            ]
            return cols[index.column()]
        elif role == MOTIF_ID_ROLE:
            return key
        elif role == Qt.BackgroundColorRole:
            return get_qt_color(index.row())
        else:
            return None

class StructureByPositionModel(QAbstractTableModel):
    def __init__(
        self,
        result : MatchMotifsResult,
        msa_to_structure_mappings : Optional[MsaToPymolStructureMap]
    ):
        super().__init__()

        # Copy the dictionary, don't trust the outside
        # forces not to modify it
        self.__result = result
        self.__keys = sorted(list(result.by_msa_position_results.keys()))
        self.__msa_to_structure_mappings = msa_to_structure_mappings

    @property
    def __front_headers(self) -> List[str]:
        return ["MSA Position", "Structure Position"]

    @property
    def __static_column_count(self) -> int:
        return len(self.__front_headers)

    def __map_to_column_to_motif(self, i : int) -> str:
        motifs = list(self.__result.motifs.keys())
        return motifs[i - len(self.__front_headers)]

    def __get_column_color(self, i : int) -> Optional[QColor]:
        if i < len(self.__front_headers):
            return None
        else:
            motif = self.__map_to_column_to_motif(i)
            return self.__result.motifs[motif].qt_color

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            headers = self.__front_headers + [f"{k} (count)" for k in self.__result.motifs.keys()]
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def rowCount(self, parent = None) -> int:
        return len(self.__result.by_msa_position_results)

    def columnCount(self, parent = None) -> int:
        return len(self.__result.motifs) + self.__static_column_count

    def __structure_position(self, position : int) -> Optional[int]:

        for msa_to_structure_mappings in viter(self.__msa_to_structure_mappings):
            for structure_pos in viter(msa_to_structure_mappings.get_pymol_structure_position(position)):
                return structure_pos

        return None

    def data(self, index: QModelIndex, role = Qt.DisplayRole):

        if not index.isValid():
            return None

        key = self.__keys[index.row()]
        entry = self.__result.by_msa_position_results[key]
        structure_position = self.__structure_position(key)

        if role == Qt.DisplayRole:
            values : List[str] = [
                str(key + 1),
                str(structure_position or ""),
            ] + [str(entry[k]) for k in self.__result.motifs.keys()]
            return values[index.column()]
        elif role == Qt.BackgroundColorRole:
            return self.__get_column_color(index.column())
        elif role == RESIDUE_ITEM_DATA_ROLE:
            return structure_position
        else:
            return None

class MatchesByNameModel(QAbstractTableModel):
    def __init__(
        self,
        result : MatchMotifsResult,
    ):
        super().__init__()

        # Copy the dictionary, don't trust the outside
        # forces not to modify it
        self.__result = result
        self.__keys = sorted(list(result.by_sequence_name_results.keys()))
        self.__motifs = list(self.__result.motifs.keys())

        self.__front_headers = ["Sequence Name"]

        self.__column_names = \
            self.__front_headers+ \
            [f"{motif} (count)" for motif in self.__motifs]

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            headers = self.__column_names
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def rowCount(self, parent = None) -> int:
        return len(self.__keys)

    def columnCount(self, parent = None) -> int:
        return len(self.__column_names)

    def __get_column_color(self, i : int) -> Optional[QColor]:
        motif_column_start = len(self.__front_headers)
        if i < motif_column_start:
            return None
        else:
            motif = self.__motifs[i - motif_column_start]
            return self.__result.motifs[motif].qt_color

    def data(self, index: QModelIndex, role = Qt.DisplayRole):

        if not index.isValid():
            return None

        key = self.__keys[index.row()]
        entry = self.__result.by_sequence_name_results[key]

        if role == Qt.DisplayRole:
            columns : List[str] = \
                [key] + \
                [str(entry.count(motif)) for motif in self.__motifs]
            return columns[index.column()]
        elif role == Qt.BackgroundColorRole:
            return self.__get_column_color(index.column())

class ColorByMotif(QWidget):

    def __init__(
        self,
        msa_context : MsaContext
    ) -> None:
        super().__init__()
        self.__ui = Ui_ColorByMotif()
        self.__ui.setupUi(self)

        self.__msa_context = msa_context

        for (name, motif_re) in KNOWN_MOTIFS.items():
            self.__ui.knownCombo.addItem(
                name,
                userData=motif_re
            )

        self.__selected_motifs : MotifsDict = {}
        self.__ui.addKnownButton.clicked.connect(self.__on_add_motif)
        self.__ui.addCustomButton.clicked.connect(self.__on_add_custom_motif)
        self.__ui.removeSelected.clicked.connect(self.__on_remove_motif)
        self.__ui.runButton.clicked.connect(self.__on_run_clicked)
        self.__ui.exportButton.clicked.connect(self.__on_export_by_sequence)

        def get_selected_structure_name():
            for selected in viter(self.__msa_context.selected_structure):
                return selected.selection
            
            return None

        self.__by_position_selection_model = pymol_residue_selection_model(
            self.__ui.structurePositionsTable,
            get_selected_structure_name,
            RESIDUE_ITEM_DATA_ROLE
        )

    def __set_motif(self, name : str, motif : re.Pattern) -> None:
        color_key = len(self.__selected_motifs)
        self.__selected_motifs[name] = PatternEntry(pattern=motif, color_index=color_key)
        self.__render_motifs_table()

    def __render_motifs_table(self):
        self.__ui.motifsTable.setModel(MotifsModel(self.__selected_motifs))
        self.__ui.motifsTable.resizeColumnsToContents()

    @pyqtSlot()
    def __on_add_motif(self):
        combo_box = self.__ui.knownCombo
        motif = combo_box.currentData(Qt.UserRole)
        self.__set_motif(combo_box.currentText(), motif)

    @pyqtSlot()
    @with_error_handler(name="__on_remove_motif")
    def __on_remove_motif(self):

        self.__remove_motifs(
            index.data(MOTIF_ID_ROLE)
            for index in self.__ui.motifsTable.selectedIndexes()
        )

    @pyqtSlot(name="__on_add_custom_motif")
    @with_error_handler()
    def __on_add_custom_motif(self):
        custom_motif_name = self.__ui.customMotifName.text().strip()
        custom_motif = self.__ui.customLineEdit.text()

        if not VALID_MOTIF_NAME.match(custom_motif_name):
            raise ValueError("The name of the motif must be only letters and digits")

        self.__add_custom_motif(custom_motif_name, custom_motif)

    def __add_custom_motif(self, name : str, pattern: str):
        self.__set_motif(name, re.compile(pattern, re.IGNORECASE))

    def __remove_motifs(self, motifs: Iterable[str]):
        for motif in motifs:
            if motif in self.__selected_motifs:
                del self.__selected_motifs[motif]

        self.__render_motifs_table()

    def __find_motifs(self, msa : Msa) -> MatchMotifsResult:

        by_position_results = {}

        def add_count(start : int, end : int, motif_name: str, count: int):

            for i in range(start, end):

                if i not in by_position_results:
                    by_position_results[i] = dict((k,0) for k in self.__selected_motifs.keys())

                by_position_results[i][motif_name] += count

        by_sequence_results = {}

        def add_sequence_match(name: str, motif : str, sequence : str, match : re.Match):

            if name not in by_sequence_results:
                by_sequence_results[name] = SequenceMatchEntry(sequence_name=name, sequence=sequence, matches=dict())

            by_sequence_results[name].add_match(motif, match)

        for name, sequence in msa.items():
            for motif_name, pattern in self.__selected_motifs.items():
                for match in pattern.finditer(sequence):
                    add_count(match.start(), match.end(), motif_name, 1)
                    add_sequence_match(name, motif_name, sequence, match)

        return MatchMotifsResult(
            by_msa_position_results=by_position_results,
            by_sequence_name_results=by_sequence_results,
            motifs=dict(self.__selected_motifs)
        )
    
    def __get_structure_mapping(self) -> Optional[MsaToPymolStructureMap]:
        for selected_structure in viter(self.__msa_context.selected_structure):
            return msa_to_pymol_structure_map(
                selected_structure,
                self.__msa_context.selected_msa_sequence_name,
                self.__msa_context.sequences,
            )


    @pyqtSlot(name="__on_run_clicked")
    @with_error_handler()
    def __on_run_clicked(self):
        results = self.__find_motifs(self.__msa_context.sequences)

        msa_to_structure_position_map = self.__get_structure_mapping()

        # Update the positions table
        self.__by_position_selection_model.setModel(
            StructureByPositionModel(
                results,
                msa_to_structure_position_map
            )
        )
        self.__ui.structurePositionsTable.resizeColumnsToContents()

        # Update the names table
        self.__ui.msaItemsTable.setModel(MatchesByNameModel(results))
        self.__ui.msaItemsTable.resizeColumnsToContents()

        self.__color_structure_by_motifs(results, msa_to_structure_position_map)

    def __color_structure_by_motifs(
        self,
        results : MatchMotifsResult,
        mappings : Optional[MsaToPymolStructureMap]
    ):

        if mappings is None:
            return

        for motif in results.motifs.keys():
            self.__color_structure_by_motif(results, motif, mappings)

    def __color_structure_by_motif(
        self,
        results : MatchMotifsResult,
        motif : str,
        msa_to_structure : MsaToPymolStructureMap
    ):
        selected_structure = self.__msa_context.selected_structure

        if selected_structure is None:
            return

        min_color_factor = 0.25
        max_color_factor = 0.99

        min_score = min(r[motif] for r in results.by_msa_position_results.values())
        max_score = max(r[motif] for r in results.by_msa_position_results.values())
        delta = max_score - min_score
        step = 0 if delta == 0 else (max_color_factor - min_color_factor) / delta
        color_range = color.color_range_scale(results.motifs[motif].rgb_color)

        for k,v in results.by_msa_position_results.items():
            resi = msa_to_structure.get_pymol_structure_position(k)
            value = v[motif]

            if resi is None or value == 0:
                continue

            scaled_score = (value * step) + min_color_factor
            pymol_color = color.to_pymol_color(color_range.get_color(scaled_score))
            pymol.cmd.color(pymol_color, f"{selected_structure.selection} and resi {resi}")

    @pyqtSlot(name="__on_export_by_sequence")
    @with_error_handler()
    def __on_export_by_sequence(self):
        filepath, _ = QFileDialog.getSaveFileName(self, "Save As", "", "CSV Files (*.csv)")

        with open(filepath, 'w') as stream:
            export_table_view_to_csv(
                self.__ui.msaItemsTable,
                stream
            )

        
