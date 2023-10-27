from Bio.Align import MultipleSeqAlignment, SeqRecord
from PyQt5.QtCore import pyqtSlot, QAbstractTableModel, QModelIndex, Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QWidget
from typing import Optional, cast, List, NamedTuple, Tuple

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

class ScoreMeta(NamedTuple):
    keep_always : bool

class SequenceScoresModel(QAbstractTableModel):

    fixed_headers = ["Always Include", "Sequence Id"]

    def __init__(
        self,
        alignment : MultipleSeqAlignment,
        scores : List[ScoreEntry]
    ):
        super().__init__()

        self.__alignment = alignment
        self.__scores = scores
        self.__score_meta = [ScoreMeta(keep_always=False) for _ in range(0, len(alignment))]

    @property
    def headers(self):
        return self.fixed_headers + [score.name for score in self.__scores]

    def rowCount(self, parent = None) -> int:
        return len(self.__alignment)

    def columnCount(self, parent = None) -> int:
        return len(self.headers)

    def set_always_keep(self, index : QModelIndex, keep : bool):
        row = index.row()
        self.__score_meta[row] = self.__score_meta[row]._replace(keep_always = keep)
        self.dataChanged.emit(index, index)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            headers = self.headers
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def __get_seq(self, index : int) -> SeqRecord:
        return cast(SeqRecord, self.__alignment[index])

    @staticmethod
    def fromat_score(value : float) -> str:
        return str(round(value, 2))

    def data(self, index: QModelIndex, role=Qt.DisplayRole):

        if not index.isValid():
            return None

        row_ix = index.row()
        col_ix = index.column()
        sequence = self.__get_seq(row_ix)
        if role == Qt.DisplayRole:
            cols : List[str] = [
                "",
                sequence.id,
                self.fromat_score(self.__scores[0].score(row_ix))
            ]
            return cols[col_ix]
        elif role == Qt.BackgroundColorRole:
            if any(not score.is_included(index.row()) for score in self.__scores):
                return QColor(255,0,0,100)
            else:
                return None
        elif role == Qt.CheckStateRole and col_ix == 0:
            return Qt.Checked
        else:
            return None

class MsaCleaner(QWidget):

    def __init__(self, ctx : Context):
        super().__init__()
        self.__ui = Ui_MsaCleaner()
        self.__ui.setupUi(self)
        self.__scores_model : Optional[SequenceScoresModel] = None

        self.__cleaners : List[Tuple[str, MsaCleanerBase]] = [
            ("Gap Divergence", ScoreByDivergence())
        ]

        for name, cleaner in self.__cleaners:
            self.__ui.cleanersWidget.addTab(cleaner, name)
            cleaner.on_score_changed.connect(self.__on_score_changed)

        self.__msa_selector = msa.msa_selector(self, self.__ui.loadMsaButton, self.__ui.selectedFileLabel)
        self.__msa_selector.msa_file_selected.connect(self.__on_msa_selected)
        self.__ui.scoresTable.clicked.connect(self.__on_sequence_score_clicked)

    @pyqtSlot()
    def __on_sequence_score_clicked(self, index : QModelIndex):

        if index.column == 0 and self.__scores_model is not None:
            self.__scores_model
        pass

    @pyqtSlot()
    def __on_score_changed(self):

        alignment = self.__msa_selector.alignment

        if alignment is None:
            return

        self.__update_table(
            SequenceScoresModel(
                alignment,
                [
                    ScoreEntry.from_result(name, result)
                    for name, cleaner in self.__cleaners
                    for result in [cleaner.score] if result is not None
                ]
            )
        )

    def __update_table(self, model : SequenceScoresModel) -> None:

        self.__scores_model = model
        self.__ui.scoresTable.setModel(model)
        self.__ui.scoresTable.resizeColumnsToContents()

    @pyqtSlot(name="__on_msa_selected")
    @with_error_handler()
    def __on_msa_selected(self):
        alignment = self.__msa_selector.alignment

        if alignment is None:
            raise Exception("The selected alignment could not be opened.")

        self.__update_table(
            SequenceScoresModel(
                alignment,
                [
                    ScoreEntry.from_result(name, result)
                    for name, cleaner in self.__cleaners
                    for result in [cleaner.score_alignment(alignment)] 
                ]
            )
        )