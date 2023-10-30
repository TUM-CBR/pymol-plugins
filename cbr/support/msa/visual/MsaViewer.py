from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QWidget
from typing import Iterable, List, Optional, Set, Tuple

from ....core import color
from ....clustal import msa
from .Ui_MsaViewer import Ui_MsaViewer

class MsaViewerModel(QAbstractTableModel):

    META_COLUMNS = ["Name"]

    def __init__(self, alignment : MultipleSeqAlignment):
        super().__init__()
        self.__alignment = alignment
        self.__residue_index = self.index_residues(alignment)
        self.__mask = []
        self.__masked_columns = []
        self.__row_mappings = list(range(0, len(alignment)))
        self.__column_mappings = list(range(0, alignment.get_alignment_length()))

    @staticmethod
    def index_residues(alignment : MultipleSeqAlignment) -> List[Set[int]]:

        return [
            set(seq_ix for seq_ix in range(0, len(alignment)) if not msa.is_blank(alignment[seq_ix][i]))
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

        self.__row_mappings = [
            next(ix for ix,_ in enumerate(self.__alignment) if ix >= row_ix and ix not in mask)
            for row_ix in range(0, self.rowCount())
        ]

        self.__column_mappings = [
            next(ix for ix in range(0, alignment_length) if ix >= col_ix and ix not in masked_columns)
            for col_ix in range(0, alignment_length - len(masked_columns))
        ]

        self.modelReset.emit()

    def __get_row_and_col(self, index : QModelIndex) -> Tuple[int, int]:
        row = self.__row_mappings[index.row()]
        column = index.column()

        meta_count = len(self.META_COLUMNS)
        if column < meta_count:
            return (row, column)

        return (row, meta_count + self.__column_mappings[column - meta_count])

    @property
    def __meta_columns_count(self) -> int:
        return len(self.META_COLUMNS)

    def __get_residue_at(self, row : int, col : int) -> str:
        return self.__alignment[row][col - len(self.META_COLUMNS)].upper()

    def __get_content_at(self, index : QModelIndex) -> str:
        row, col = self.__get_row_and_col(index)
        if col == 0:
            return self.__alignment[row].id
        else:
            return self.__get_residue_at(row, col)

    def __get_color_at(self, index : QModelIndex) -> Optional[QColor]:

        if index.column() < self.__meta_columns_count:
            return None

        row, col = self.__get_row_and_col(index)
        resi = self.__get_residue_at(row, col)
        m_rgb = color.residue_letter_colors.get(resi)

        if m_rgb is None:
            return None

        return QColor(*m_rgb)

    def rowCount(self, parent = None) -> int:
        return len(self.__alignment) - len(self.__mask)

    def columnCount(self, parent = None) -> int:
        return len(self.META_COLUMNS) + self.__alignment.get_alignment_length() - len(self.__masked_columns)

    def data(self, index: QModelIndex, role=Qt.DisplayRole):

        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            return self.__get_content_at(index)
        if role == Qt.BackgroundColorRole:
            return self.__get_color_at(index)
        else:
            return None

class MsaViewer(QWidget):

    def __init__(self):
        super().__init__()
        self.__ui = Ui_MsaViewer()
        self.__ui.setupUi(self)
        self.__model : Optional[MsaViewerModel] = None

    def __adjust_columns_size(self) -> None:
        model = self.__model

        if model is None:
            return

        self.__ui.msaTable.resizeColumnToContents(0)
        for i in range(1, model.columnCount()):
            self.__ui.msaTable.setColumnWidth(i, 10)

    def mask_sequences(self, rows : Iterable[int]) -> None:

        if self.__model is None:
            return

        self.__model.mask_sequences(rows)
        self.__adjust_columns_size()

    def set_alignment(self, alignment : MultipleSeqAlignment) -> None:

        self.__model = model = MsaViewerModel(alignment)
        self.__ui.msaTable.setModel(model)
        self.__adjust_columns_size()

