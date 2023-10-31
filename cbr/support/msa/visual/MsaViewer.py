from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from enum import Enum
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt, pyqtSlot
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QWidget
from typing import Iterable, List, Optional, Set, Tuple

from ....core import color
from ....clustal import msa
from .Ui_MsaViewer import Ui_MsaViewer

class MaskPositionMode(Enum):
    HIGHLIGHT = 0
    HIDE = 1

class MsaViewerModel(QAbstractTableModel):

    META_COLUMNS = ["Name"]
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

        self.__masked_columns = masked_columns = [
            column
            for column in range(0, alignment_length)
            for residue_index in [self.__residue_index[column]] if len(residue_index.difference(mask)) == 0
        ]

        self.__masked_row_mappings = [
            row_ix + sum(1 for ix in self.__mask if row_ix >= ix)
            for row_ix in range(0, self.masked_alignment_seqs)
        ]

        self.__masked_column_mappings = [
            col_ix + sum(1 for ix in masked_columns if col_ix >= ix)
            for col_ix in range(0, self.masked_alignment_length)
        ]

        self.modelReset.emit()

    def __get_row_mapping(self, row : int) -> int:
        if self.__mask_position_mode == MaskPositionMode.HIDE:
            return self.__masked_row_mappings[row]

        return row

    def __get_column_mapping(self, seq_column_ix : int) -> int:

        if self.__mask_position_mode == MaskPositionMode.HIDE:
            return self.__masked_column_mappings[seq_column_ix]

        return seq_column_ix

    def __get_row_and_col(self, index : QModelIndex) -> Tuple[int, int]:
        row = self.__get_row_mapping(index.row())
        column = index.column()

        if self.__is_meta(index):
            return (row, column)

        row, column = self.__get_structure_position(index)
        return (row, len(self.META_COLUMNS) + column)

    def __get_structure_position(self, index : QModelIndex) -> Tuple[int, int]:
        assert not self.__is_meta(index) >= len(self.META_COLUMNS), "The index is a metadata column"

        return (index.row(), self.__get_column_mapping(index.column() - len(self.META_COLUMNS)))

    @property
    def __meta_columns_count(self) -> int:
        return len(self.META_COLUMNS)

    def __get_residue_at(self, row : int, col : int) -> str:
        return self.__alignment[row][col - len(self.META_COLUMNS)].upper() # type: ignore[reportGeneralTypeIssues]

    def __get_content_at(self, index : QModelIndex) -> str:
        row, col = self.__get_row_and_col(index)
        if col == 0:
            return self.__alignment[row].id # type: ignore[reportGeneralTypeIssues]
        else:
            return self.__get_residue_at(row, col)

    def __is_masked(self, index : QModelIndex) -> bool:

        if self.__is_meta(index):
            return False

        row, col = self.__get_structure_position(index)
        return row in self.__mask or col in self.__masked_columns

    def __get_color_at(self, index : QModelIndex) -> Optional[QColor]:

        if self.__is_meta(index):
            return None

        if self.__mask_position_mode == MaskPositionMode.HIGHLIGHT and self.__is_masked(index):
            return self.BLACK

        return self.__get_resiude_color(index)

    RESIDUE_COLORS = dict((resi, QColor(*rgb)) for resi, rgb in color.residue_letter_colors.items())

    def __get_resiude_color(self, index: QModelIndex) -> Optional[QColor]:
        row, col = self.__get_row_and_col(index)
        resi = self.__get_residue_at(row, col)

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
        self.__ui.msaTable.setModel(model)
        self.__adjust_columns_size()

    def get_new_alignment(self) -> MultipleSeqAlignment:
        """Get the resulting alignment when the selected rows are masked"""

        if self.__model is None:
            raise ValueError("No alignment has been provided.")

        return self.__model.get_masked_alignment()

