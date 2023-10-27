from Bio.Align import MultipleSeqAlignment, SeqRecord
from PyQt5.QtCore import pyqtSlot, QAbstractTableModel, QModelIndex, Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QWidget
from typing import cast, List, NamedTuple, Tuple

from ...core.Context import (Context)
from ...core.Qt.QtWidgets import with_error_handler
from ...support import msa

from .Ui_MsaCleaner import Ui_MsaCleaner
from .MsaCleanerResult import MsaCleanerBase, MsaCleanerResult
from .ScoreByDivergence import ScoreByDivergence

class ScoreEntry(NamedTuple):
    name : str
    scores : List[float]
    treshold : Tuple[float, float]

    def is_included(self, i :int) -> bool:
        low, high = self.treshold
        score = self.scores[i]
        return score >= low and score <= high

    @staticmethod
    def from_result(name : str, result : MsaCleanerResult) -> 'ScoreEntry':
        return ScoreEntry(
            name,
            result.scores,
            result.treshold
        )

    def score(self, i : int) -> float:
        return self.scores[i]

class SequenceScoresModel(QAbstractTableModel):

    fixed_headers = ["Sequence Id"]

    def __init__(
        self,
        alignment : MultipleSeqAlignment,
        scores : List[ScoreEntry]
    ):
        super().__init__()

        self.__alignment = alignment
        self.__scores = scores

    @property
    def headers(self):
        return self.fixed_headers + [score.name for score in self.__scores]

    def rowCount(self, parent = None) -> int:
        return len(self.__alignment)

    def columnCount(self, parent = None) -> int:
        return len(self.headers)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            headers = self.headers
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def __get_seq(self, index : int) -> SeqRecord:
        return cast(SeqRecord, self.__alignment[index])

    def data(self, index: QModelIndex, role=Qt.DisplayRole):

        if not index.isValid():
            return None

        row_ix = index.row()
        col_ix = index.column()
        sequence = self.__get_seq(row_ix)
        if role == Qt.DisplayRole:
            cols : List[str] = [
                sequence.id,
                str(self.__scores[0].score(row_ix))
            ]
            return cols[col_ix]
        elif role == Qt.BackgroundColorRole:
            if any(not score.is_included(index.row()) for score in self.__scores):
                return QColor(255,0,0,100)
            else:
                return None
        else:
            return None

class MsaCleaner(QWidget):

    def __init__(self, ctx : Context):
        super().__init__()
        self.__ui = Ui_MsaCleaner()
        self.__ui.setupUi(self)

        self.__cleaners : List[Tuple[str, MsaCleanerBase]] = [
            ("Gap Divergence", ScoreByDivergence())
        ]

        for name, cleaner in self.__cleaners:
            self.__ui.cleanersWidget.addTab(cleaner, name)

        self.__msa_selector = msa.msa_selector(self, self.__ui.loadMsaButton, self.__ui.selectedFileLabel)
        self.__msa_selector.msa_file_selected.connect(self.__on_msa_selected)

    @pyqtSlot(name="__on_msa_selected")
    @with_error_handler()
    def __on_msa_selected(self):
        alignment = self.__msa_selector.alignment

        if alignment is None:
            raise Exception("The selected alignment could not be opened.")

        self.__ui.scoresTable.setModel(
            SequenceScoresModel(
                alignment,
                [
                    ScoreEntry.from_result(name, result)
                    for name, cleaner in self.__cleaners
                    for result in [cleaner.score_alignment(alignment)] 
                ]
            )
        )