from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from enum import Enum
import pymol
from PyQt5.QtCore import QAbstractTableModel, QItemSelection, QItemSelectionModel, QModelIndex, Qt, pyqtSlot
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QDialog, QWidget
from typing import Dict, Iterable, List, NamedTuple, Optional, Set

from ....core import color
from ....core.pymol.structure import StructureSelection
from ....clustal import msa
from ....control import viter, update_key
from .MsaStructureSelector import MsaStructureSelector
from ..structure import msa_to_pymol_structure_map, MsaToPymolStructureMap
from .Ui_MsaViewer import Ui_MsaViewer

class MaskPositionMode(Enum):
    HIGHLIGHT = 0
    HIDE = 1

def mask_offset(
    mask : Set[int],
    length : int
):
    masked_mappings = []
    col_ix = 0

    while(col_ix < length):
        ix = 0 if len(masked_mappings) == 0 else 1 + masked_mappings[-1]
        while(ix in mask):
            ix += 1

        assert ix not in mask

        masked_mappings.append(ix)
        col_ix += 1

    return masked_mappings

class SequenceMeta(NamedTuple):
    sequence : SeqRecord
    structure : Optional[StructureSelection]
    structure_to_sequence : Optional[MsaToPymolStructureMap]

    @staticmethod
    def create(seq : SeqRecord):
        return SequenceMeta(seq, None, None)

    def get_structure_position(self, pos: int) -> Optional[int]:
        if self.structure_to_sequence is None:
            return None
        else:
            self.structure_to_sequence.get_pymol_structure_position(pos)

    def update_structure(
        self,
        msa : MultipleSeqAlignment,
        structure: Optional[StructureSelection]
    ):

        if structure is None:
            return self._replace(structure = None, structure_to_sequence = None)

        mapping = msa_to_pymol_structure_map(
            structure,
            self.sequence.id,
            msa,
        )

        return self._replace(
            structure = structure,
            structure_to_sequence = mapping
        )

class MsaViewerModel(QAbstractTableModel):

    META_COLUMNS = ["Name", "Structure"]
    STRUCTURE_COLUMN = META_COLUMNS.index("Structure")
    BLACK = QColor(0,0,0)

    def __init__(self, alignment : MultipleSeqAlignment):
        super().__init__()
        self.__alignment = alignment
        self.__residue_index = self.index_residues(alignment)
        self.__mask = set([])
        self.__masked_columns = []
        self.__masked_row_mappings = list(range(0, len(alignment)))
        self.__masked_column_mappings = list(range(0, alignment.get_alignment_length()))
        self.__mask_position_mode = MaskPositionMode.HIDE
        self.__sequences_meta = [SequenceMeta.create(seq) for seq in alignment]

    def get_positions_from_selection(self, selection : QItemSelection) -> Dict[str, List[int]]:

        result = dict()
        for index in selection.indexes():
            if self.__is_meta(index):
                continue

            msa_row = self.__get_row_mapping(index.row())
            msa_col = self.__get_column_mapping(index.column())
            meta = self.__get_sequence_meta(msa_row)

            if meta.structure is None:
                continue

            sele_str = meta.structure.selection

            with update_key(result, sele_str, []) as update:
                update.value = update.value + [meta.get_structure_position(msa_col)]

        return result

    def __get_sequence_meta(self, col : int) -> SequenceMeta:
        return self.__sequences_meta[col]

    def set_structure(self, index: QModelIndex, structure: Optional[StructureSelection]):
        assert index.column() == self.STRUCTURE_COLUMN
        row = index.row()
        self.__sequences_meta[row] = self.__sequences_meta[row].update_structure(
            self.__alignment,
            structure = structure
        )
        self.dataChanged.emit(index, index, [Qt.DisplayRole])

    def __get_headers(self):

        return self.META_COLUMNS \
            + [str(i + 1) for i in range(0, self.__alignment.get_alignment_length())]

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            headers = self.__get_headers()
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def get_masked_alignment(self) -> MultipleSeqAlignment:

        return MultipleSeqAlignment(
            SeqRecord(
                Seq("".join(c for i,c in enumerate(seq) if i not in self.__masked_columns)),
                id = seq.id
            )
            for row_ix, seq in enumerate(self.__alignment) if row_ix not in self.__mask
        )

    @staticmethod
    def index_residues(alignment : MultipleSeqAlignment) -> List[Set[int]]:

        return [
            set(seq_ix for seq_ix,seq in enumerate(alignment) if not msa.is_blank(seq[i])) # type: ignore[reportGeneralTypeIssues]
            for i in range(0, alignment.get_alignment_length())
        ]

    def mask_sequences(self, mask : Iterable[int]) -> None:
        """Exclude the sequences at the given indexes from display"""

        self.__mask = mask = set(mask)
        alignment_length = self.__alignment.get_alignment_length()

        self.__masked_columns = masked_columns = set(
            column
            for column in range(0, alignment_length)
            for residue_index in [self.__residue_index[column]] if len(residue_index.difference(mask)) == 0
        )

        self.__masked_row_mappings = mask_offset(mask, self.masked_alignment_seqs)
        self.__masked_column_mappings = mask_offset(masked_columns, self.masked_alignment_length)

        self.modelReset.emit()

    def __get_row_mapping(self, row : int) -> int:
        if self.__mask_position_mode == MaskPositionMode.HIDE:
            return self.__masked_row_mappings[row]

        return row

    def __get_column_mapping(self, seq_column_ix : int) -> int:

        seq_column_ix = seq_column_ix - self.__meta_columns_count
        if self.__mask_position_mode == MaskPositionMode.HIDE:
            return self.__masked_column_mappings[seq_column_ix]

        return seq_column_ix

    @property
    def __meta_columns_count(self) -> int:
        return len(self.META_COLUMNS)

    def __get_residue_at(self, row : int, col : int) -> str:
        return self.__alignment[row][col].upper() # type: ignore[reportGeneralTypeIssues]

    def __get_meta_content(self, row: int, col : int) -> str:
        structure = self.__sequences_meta[row].structure
        content = [
            self.__alignment[row].id, # type: ignore[reportGeneralTypeIssues],
            structure.show() if structure else "<double click to select>"
        ]

        return content[col]

    def __get_content_at(self, index : QModelIndex) -> str:
        row = self.__get_row_mapping(index.row())
        if self.__is_meta(index):
            return self.__get_meta_content(row, index.column())

        col = self.__get_column_mapping(index.column())
        return self.__get_residue_at(row, col)

    def __is_masked(self, index : QModelIndex) -> bool:

        if self.__is_meta(index):
            return False

        row = index.row()
        col = index.column() - self.__meta_columns_count

        return row in self.__mask or col in self.__masked_columns

    def __get_color_at(self, index : QModelIndex) -> Optional[QColor]:

        if self.__is_meta(index):
            return None

        if self.__mask_position_mode == MaskPositionMode.HIGHLIGHT and self.__is_masked(index):
            return self.BLACK

        return self.__get_resiude_color(index)

    RESIDUE_COLORS = dict((resi, QColor(*rgb)) for resi, rgb in color.residue_letter_colors.items())

    def __get_resiude_color(self, index: QModelIndex) -> Optional[QColor]:

        assert not self.__is_meta(index), "The meta columns don't have a residue color"

        resi = self.__get_content_at(index)

        return self.RESIDUE_COLORS.get(resi)

    def rowCount(self, parent = None) -> int:

        if self.__mask_position_mode == MaskPositionMode.HIDE:
            return len(self.__alignment) - len(self.__mask)
        else:
            return len(self.__alignment)

    def columnCount(self, parent = None) -> int:
        if self.__mask_position_mode == MaskPositionMode.HIDE:
            return len(self.META_COLUMNS) + self.masked_alignment_length
        else:
            return len(self.META_COLUMNS) + self.__alignment.get_alignment_length()

    @property
    def masked_alignment_length(self) -> int:
        return self.__alignment.get_alignment_length() - len(self.__masked_columns)

    @property
    def masked_alignment_seqs(self) -> int:
        return len(self.__alignment) - len(self.__mask)

    def set_mask_position_mode(self, mode: MaskPositionMode):
        self.__mask_position_mode = mode
        self.modelReset.emit()

    def __is_meta(self, index: QModelIndex) -> bool:
        return index.column() < len(self.META_COLUMNS)

    WHITE = QColor(255, 255, 255)

    def __get_text_color_at(self, index: QModelIndex) -> Optional[QColor]:

        if not self.__is_meta(index) \
            and self.__mask_position_mode == MaskPositionMode.HIGHLIGHT \
            and self.__is_masked(index):
            return self.__get_resiude_color(index) or self.WHITE

        return None

    def data(self, index: QModelIndex, role=Qt.DisplayRole):

        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            return self.__get_content_at(index)
        if role == Qt.BackgroundColorRole:
            return self.__get_color_at(index)
        if role == Qt.TextColorRole:
            return self.__get_text_color_at(index)
        else:
            return None

MASK_POSITION_MODE = {
    "Hide": MaskPositionMode.HIDE,
    "Highlight": MaskPositionMode.HIGHLIGHT
}

class MsaViewer(QWidget):

    def __init__(self):
        super().__init__()
        self.__ui = Ui_MsaViewer()
        self.__ui.setupUi(self)
        self.__model : Optional[MsaViewerModel] = None

        self.__ui.maskCombo.addItems(MASK_POSITION_MODE.keys())
        self.__ui.maskCombo.currentIndexChanged.connect(self.__on_mask_combo_changed)
        self.__ui.msaTable.doubleClicked.connect(self.__select_structure)
        self.__select_structure_dialog = MsaStructureSelector(self)
        self.__ui.msaTable.setSelectionModel(QItemSelectionModel())

    @pyqtSlot(QItemSelection, QItemSelection)
    def __on_selection(self, selected: QItemSelection, deselected: QItemSelection):

        for model in viter(self.__model):
            for model_sele, resi in model.get_positions_from_selection(selected):
                resi_sele = " or ".join(f"resi {i}" for i in resi)
                pymol.cmd.select(
                    "sele (msa cleaner)",
                    f"{model_sele} and {resi_sele}"
                )

    @pyqtSlot()
    def __on_mask_combo_changed(self):

        if self.__model is None:
            return

        self.__model.set_mask_position_mode(MASK_POSITION_MODE[self.__ui.maskCombo.currentText()])
        self.__adjust_columns_size()

    def __adjust_columns_size(self) -> None:
        model = self.__model

        if model is None:
            return

        self.__ui.msaTable.resizeColumnToContents(0)
        for i in range(1, model.columnCount()):
            self.__ui.msaTable.setColumnWidth(i, 10)

    @pyqtSlot(QModelIndex)
    def __select_structure(self, index : QModelIndex):

        for model in viter(self.__model):
            if (self.__select_structure_dialog.exec() == QDialog.Accepted):
                model.set_structure(
                    index,
                    self.__select_structure_dialog.current_selection()
                )

    @property
    def masked_alignment_length(self) -> int:

        if self.__model is None:
            return 0

        return self.__model.masked_alignment_length

    @property
    def masked_alignment_seqs(self) -> int:

        if self.__model is None:
            return 0

        return self.__model.masked_alignment_seqs

    def mask_sequences(self, rows : Iterable[int]) -> None:

        if self.__model is None:
            return

        self.__model.mask_sequences(rows)
        self.__adjust_columns_size()

    def set_alignment(self, alignment : MultipleSeqAlignment) -> None:

        self.__model = model = MsaViewerModel(alignment)
        selection_model = QItemSelectionModel(model)
        selection_model.selectionChanged.connect(self.__on_selection)
        self.__ui.msaTable.setModel(model)
        self.__ui.msaTable.setSelectionModel(selection_model)
        self.__adjust_columns_size()

    def get_new_alignment(self) -> MultipleSeqAlignment:
        """Get the resulting alignment when the selected rows are masked"""

        if self.__model is None:
            raise ValueError("No alignment has been provided.")

        return self.__model.get_masked_alignment()

