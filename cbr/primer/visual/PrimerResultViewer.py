import math
from PyQt5.QtWidgets import QWidget
import textwrap
from typing import NamedTuple

from ..Operations import PrimerResult
from .Ui_PrimerResultViewer import Ui_PrimerResultViewer

class PrimerResultContext(NamedTuple):
    result : PrimerResult
    original_sequence : str

class PrimerResultViewer(QWidget):

    def __init__(self):

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

        i_left = original_sequence.find(primer_result.left_primer)
        i_right = original_sequence.find(primer_result.right_primer)

        if i_left + len(primer_result.inner_seq) != i_right:
            raise Exception("The provided sequence does not match the provided primers")

        for i,base in enumerate(original_sequence):
            pass

