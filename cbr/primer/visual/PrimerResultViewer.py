import math
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QTableWidgetItem, QWidget
from typing import NamedTuple

from ..Operations import PrimerResult
from .Ui_PrimerResultViewer import Ui_PrimerResultViewer

class PrimerResultContext(NamedTuple):
    result : PrimerResult
    original_sequence : str

PRIMER_COLOR = QColor(0,0,100,100)
LEFT_PRIMER_COLOR = PRIMER_COLOR
RIGHT_PRIMER_COLOR = PRIMER_COLOR
CODON_COLOR = QColor(0,100,0,100)

class PrimerResultViewer(QWidget):

    def __init__(self, father : QWidget):

        super(PrimerResultViewer, self).__init__(father)

        self.__ui = Ui_PrimerResultViewer()
        self.__ui.setupUi(self)

    @property
    def __line_size(self) -> int:
        return 60

    def set_result(self, result : PrimerResultContext):
        table = self.__ui.resultTable
        table.clearContents()

        original_sequence = result.original_sequence
        primer_result = result.result

        num_sections = math.ceil(len(original_sequence) / self.__line_size)
        row_count = num_sections * 3
        col_count = self.__line_size
        table.setRowCount(row_count)
        table.setColumnCount(col_count)

        i_left = original_sequence.find(primer_result.left_primer)
        i_left_end = i_left + len(primer_result.left_primer)
        c_left = primer_result.c_left_primer
        i_right = original_sequence.find(primer_result.right_primer)
        i_right_end = i_right + len(primer_result.right_primer)
        c_right = primer_result.c_right_primer
        c_codon = primer_result.c_inner

        if i_left_end + len(primer_result.inner_seq) != i_right:
            raise Exception("The provided sequence does not match the provided primers")

        for i,base in enumerate(original_sequence):
            row = math.floor(i / col_count) * 3
            col = i % col_count
            seq_item = QTableWidgetItem(base)

            if i >= i_left and i < i_left_end:
                seq_item.setBackground(LEFT_PRIMER_COLOR)
                l_primer_item = QTableWidgetItem(c_left[i - i_left])
                l_primer_item.setBackground(LEFT_PRIMER_COLOR)
                table.setItem(row + 1, col, l_primer_item)

            if i >= i_right and i < i_right_end:
                seq_item.setBackground(RIGHT_PRIMER_COLOR)
                r_primer_item = QTableWidgetItem(c_right[i - i_right])
                r_primer_item.setBackground(RIGHT_PRIMER_COLOR)
                table.setItem(row + 1, col, r_primer_item)

            if i >= i_left_end and i < i_right:
                seq_item.setBackground(CODON_COLOR)
                codon_item = QTableWidgetItem(c_codon[i - i_left_end])
                codon_item.setBackground(CODON_COLOR)
                table.setItem(row + 1, col, codon_item)

            table.setItem(row, col, seq_item)

        table.resizeColumnsToContents()

