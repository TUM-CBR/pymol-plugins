"""This scoring function is designed to eliminate sequences that create long inserts in the MSA.
It has two parameters. The first one "Continuity Treshold" determines the percentage of sequences
that must have gaps at a given position for two consecutive residues to be considered an "insert".
The second parameter "Score Treshold" controls how many sequences should be dropped based on the
score the received."""

from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import pyqtSlot
from typing import Optional, Tuple

from ...core.Qt.QtWidgets import throttle
from ..cleanup import score_inserts
from .MsaCleanerResult import MsaCleanerBase, MsaCleanerResult
from .Ui_ScoreByLongInserts import Ui_ScoreByLongInserts

class ScoreByLongInserts(MsaCleanerBase):

    def __init__(self):
        super().__init__()
        self.__ui = Ui_ScoreByLongInserts()
        self.__ui.setupUi(self)
        self.__score : Optional[MsaCleanerResult] = None

        assert(__doc__)
        self.__ui.descriptionLabel.setText(__doc__.replace("\n", " "))
        self.__alignment : Optional[MultipleSeqAlignment] = None

        self.__ui.continuityTresholdSlider.valueChanged.connect(self.__update_score)
        self.__ui.scoreTresholdSlider.valueChanged.connect(self.__update_treshold)

    @pyqtSlot(name="__update_score")
    @throttle(1000)
    def __update_score(self):

        if self.__alignment is None:
            return

        self.__score = MsaCleanerResult(
            scores = score_inserts(
                self.__alignment,
                self.__continuity_treshold
            ),
            treshold=self.__treshold
        )

        self.on_score_changed.emit()

    @pyqtSlot(name="__update_score")
    @throttle(1000)
    def __update_treshold(self):
        if self.__score is None:
            return

        self.__score = self.__score._replace(treshold = self.__treshold)
        self.on_score_changed.emit()

    @property
    def __treshold(self) -> Tuple[float, float]:
        return (0, self.__ui.scoreTresholdSlider.value() / 100)

    @property
    def __continuity_treshold(self) -> float:
        return self.__ui.continuityTresholdSlider.value() / 100

    def score_alignment(self, alignment : MultipleSeqAlignment) -> MsaCleanerResult:
        self.__alignment = alignment
        msa_score = score_inserts(alignment, self.__continuity_treshold)
        result = MsaCleanerResult(
            scores=msa_score,
            treshold=self.__treshold
        )

        self.__score = result
        return result

    @property
    def score(self) -> Optional[MsaCleanerResult]:
        return self.__score