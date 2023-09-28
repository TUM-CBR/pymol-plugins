import math
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QTableWidgetItem, QWidget
from typing import Callable, NamedTuple, Optional, TypeVar

from ..data import DesignPrimersResults, PrimerResult
from .Ui_PrimerResultViewer import Ui_PrimerResultViewer

class PrimerResultContext(NamedTuple):
    result : Optional[PrimerResult]
    source : DesignPrimersResults

    @property
    def original_sequence(self) -> str:
        return self.source.plasmid

T = TypeVar('T')

class PrimerRegionMeta(NamedTuple):

    class Fields(NamedTuple):
        i_left : int  
        i_left_end : int
        c_left : str
        i_right : int
        i_right_end : int
        c_right : str
        c_codon : str

        def left_with_offset(self, i : int) -> Optional[str]:
            if self.i_left <= i and i < self.i_left_end:
                return self.c_left[i - self.i_left]
            
            return None

        def right_with_offset(self, i : int) -> Optional[str]:
            if self.i_right <= i and i < self.i_right_end:
                return self.c_right[i - self.i_right]

            return None

        def inner_with_offset(self, i : int) -> Optional[str]:
            if self.i_left_end <= i and i < self.i_right:
                return self.c_codon[i - self.i_left_end]

            return None

    fields : Optional[Fields]

    @staticmethod
    def from_result(result : PrimerResultContext) -> 'PrimerRegionMeta':

        primer_result = result.result
        if primer_result is None:
            return PrimerRegionMeta(fields=None)

        original_sequence = result.original_sequence

        i_left = original_sequence.find(primer_result.left_primer)
        i_left_end = i_left + len(primer_result.left_primer)
        c_left = primer_result.c_left_primer
        i_right = original_sequence.find(primer_result.right_primer)
        i_right_end = i_right + len(primer_result.right_primer)

        if i_left_end + len(primer_result.inner_seq) != i_right:
            raise Exception("The provided sequence does not match the provided primers")

        return PrimerRegionMeta(
            fields = PrimerRegionMeta.Fields (
                i_left = i_left,
                i_left_end = i_left_end,
                c_left = c_left,
                i_right = i_right,
                i_right_end = i_right_end,
                c_right = primer_result.c_right_primer,
                c_codon = primer_result.c_inner
            )
        )

    def __with_fields(
        self,
        check : 'Callable[[PrimerRegionMeta.Fields], T]',
        default : T
    ):
        if self.fields is None:
            return default

        return check(self.fields)

    def left_with_offset(self, i : int) -> Optional[str]:
        return self.__with_fields(
            lambda fields: fields.left_with_offset(i),
            None
        )

    def right_with_offset(self, i : int) -> Optional[str]:
        return self.__with_fields(
            lambda fields: fields.right_with_offset(i),
            None
        )

    def inner_with_offset(self, i : int) -> Optional[str]:
        return self.__with_fields(
            lambda fields: fields.inner_with_offset(i),
            None
        )

PRIMER_COLOR = QColor(0,0,100,100)
LEFT_PRIMER_COLOR = PRIMER_COLOR
RIGHT_PRIMER_COLOR = PRIMER_COLOR
CODON_COLOR = QColor(0,100,0,100)
ROWS_MULTIPLIER = 3

class PrimerResultViewer(QWidget):

    def __init__(self, father : QWidget):

        super(PrimerResultViewer, self).__init__(father)

        self.__ui = Ui_PrimerResultViewer()
        self.__ui.setupUi(self)

    @property
    def __line_size(self) -> int:
        return 60

    @property
    def __rows_per_line(self) -> int:
        return ROWS_MULTIPLIER

    def set_result(self, result : PrimerResultContext):
        table = self.__ui.resultTable
        table.clearContents()

        original_sequence = result.original_sequence
        primer_result = result.result

        num_sections = math.ceil(len(original_sequence) / self.__line_size)
        row_count = num_sections * 3

    
        col_count = self.__line_size
        rows_per_line = self.__rows_per_line
        table.setRowCount(row_count)
        table.setColumnCount(col_count)

        region_meta = PrimerRegionMeta.from_result(result)

        for i,base in enumerate(original_sequence):
            row = math.floor(i / col_count) * rows_per_line
            col = i % col_count
            seq_item = QTableWidgetItem(base)

            primer_base = region_meta.left_with_offset(i)
            if primer_base is not None:
                seq_item.setBackground(LEFT_PRIMER_COLOR)
                l_primer_item = QTableWidgetItem(primer_base)
                l_primer_item.setBackground(LEFT_PRIMER_COLOR)
                table.setItem(row + 1, col, l_primer_item)

            primer_base = region_meta.right_with_offset(i)
            if primer_base is not None:
                seq_item.setBackground(RIGHT_PRIMER_COLOR)
                r_primer_item = QTableWidgetItem(primer_base)
                r_primer_item.setBackground(RIGHT_PRIMER_COLOR)
                table.setItem(row + 1, col, r_primer_item)

            primer_base = region_meta.inner_with_offset(i)
            if primer_base is not None:
                seq_item.setBackground(CODON_COLOR)
                codon_item = QTableWidgetItem(primer_base)
                codon_item.setBackground(CODON_COLOR)
                table.setItem(row + 1, col, codon_item)

            table.setItem(row, col, seq_item)

        table.resizeColumnsToContents()

