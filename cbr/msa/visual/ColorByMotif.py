import pymol
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, pyqtSlot, Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QWidget
import re
from typing import Dict, Iterable, Iterator, List, NamedTuple, Optional

from ...core import color
from ...core.color import RgbColor
from ...core.Qt.QtWidgets import with_error_handler
from ...support.msa import Msa

from .MsaContext import MsaContext
from .support import msa_to_structure_position_map, sequence_to_structure_position_map
from .Ui_ColorByMotif import Ui_ColorByMotif

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

class PatternEntry(NamedTuple):
    pattern : re.Pattern
    color_index : int

    @property
    def rgb_color(self) -> RgbColor:
        return color.get_color_rgb(self.color_index)

    @property
    def qt_color(self) -> QColor:
        return QColor(*color.get_qt_color(self.color_index))

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
            return QColor(*color.get_qt_color(index.row()))
        else:
            return None

class StructureByPositionModel(QAbstractTableModel):
    def __init__(
        self,
        result : MatchMotifsResult
    ):
        super().__init__()

        # Copy the dictionary, don't trust the outside
        # forces not to modify it
        self.__result = result
        self.__keys = sorted(list(result.by_msa_position_results.keys()))

    @property
    def __front_headers(self) -> List[str]:
        return ["MSA Index"]

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

    def data(self, index: QModelIndex, role = Qt.DisplayRole):

        if not index.isValid():
            return None

        key = self.__keys[index.row()]
        entry = self.__result.by_msa_position_results[key]

        if role == Qt.DisplayRole:
            values : List[str] = [
                str(key + 1)
            ] + [str(entry[k]) for k in self.__result.motifs.keys()]
            return values[index.column()]
        elif role == Qt.BackgroundColorRole:
            return self.__get_column_color(index.column())
        else:
            return None

class MatchesByNameModel(QAbstractTableModel):
    def __init__(
        self,
        result : MatchMotifsResult
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
        self.__set_motif(name, re.compile(pattern))

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

    @pyqtSlot(name="__on_run_clicked")
    @with_error_handler()
    def __on_run_clicked(self):
        results = self.__find_motifs(self.__msa_context.sequences)

        # Update the positions table
        self.__ui.structurePositionsTable.setModel(StructureByPositionModel(results))
        self.__ui.structurePositionsTable.resizeColumnsToContents()

        # Update the names table
        self.__ui.msaItemsTable.setModel(MatchesByNameModel(results))
        self.__ui.msaItemsTable.resizeColumnsToContents()

        self.__color_structure_by_motifs(results)

    def __color_structure_by_motifs(self, results : MatchMotifsResult):

        if self.__msa_context.selected_structure is None:
            return

        mappings = msa_to_structure_position_map(
            self.__msa_context.selected_msa_sequence_name,
            self.__msa_context.sequences,
            self.__msa_context.get_structure_sequence()
        )

        for motif in results.motifs.keys():
            self.__color_structure_by_motif(results, motif, mappings)

    def __color_structure_by_motif(
        self,
        results : MatchMotifsResult,
        motif : str,
        msa_to_structure : List[Optional[int]]
    ):
        selected_structure = self.__msa_context.selected_structure

        if selected_structure is None:
            return

        min_color_factor = 0.5
        max_color_factor = 0.99

        min_score = min(r[motif] for r in results.by_msa_position_results.values())
        max_score = max(r[motif] for r in results.by_msa_position_results.values())
        step = (max_color_factor - min_color_factor) / (max_score - min_score)
        structure_mappings = sequence_to_structure_position_map(selected_structure)
        color_range = color.color_range_scale(results.motifs[motif].rgb_color)

        for k,v in results.by_msa_position_results.items():
            structure_pos = msa_to_structure[k]

            if structure_pos is None:
                continue

            resi = structure_mappings[structure_pos]
            scaled_score = (v[motif] * step) + min_color_factor
            pymol_color = color.to_pymol_color(color_range.get_color(scaled_score))
            pymol.cmd.color(pymol_color, f"{selected_structure.selection} and resi {resi}")