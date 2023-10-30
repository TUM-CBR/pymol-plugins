"""This score is defined as the number of gaps introduced by a sequence devided by the total number of
of gaps in the multiple sequence alignment. The higher the score means that a greater percentage of the
gaps are introduced by the sequence.
"""

from Bio.Align import MultipleSeqAlignment
from typing import Optional, Tuple

from PyQt5.QtCore import pyqtSlot

from ..cleanup import score_by_gap_divergence
from .Ui_ScoreByDivergence import Ui_ScoreByDivergence
from .MsaCleanerResult import MsaCleanerBase, MsaCleanerResult

class ScoreByDivergence(MsaCleanerBase):

    def __init__(self):
        super().__init__()
        self.__ui = Ui_ScoreByDivergence()
        self.__ui.setupUi(self)
        assert(__doc__)
        self.__ui.descriptionLabel.setText(__doc__.replace("\n", " "))
        self.__score : Optional[MsaCleanerResult] = None
        self.__ui.divergenceSlider.valueChanged.connect(self.__on_range_changed)
        self.__update_score_label()

    def __update_score_label(self):
        self.__ui.valueLabel.setText(str(self.__treshold[1]))

    @pyqtSlot()
    def __on_range_changed(self):

        if self.__score is None:
            return

        self.__score = self.__score._replace(treshold = self.__treshold)
        self.__update_score_label()
        self.on_score_changed.emit()

    @property
    def score(self) -> Optional[MsaCleanerResult]:
        return self.__score

    @property
    def __treshold(self) -> Tuple[float, float]:
        return (0, self.__ui.divergenceSlider.value() / 100)

    def score_alignment(self, alignment : MultipleSeqAlignment) -> MsaCleanerResult:
        msa_score = score_by_gap_divergence(alignment)
        result = MsaCleanerResult(
            scores=msa_score,
            treshold=self.__treshold
        )

        self.__score = result
        return result
