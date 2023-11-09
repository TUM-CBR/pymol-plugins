"""This scoring function is designed to eliminate sequences that create long inserts in the MSA.
It has two parameters. The first one "Continuity Treshold" determines the percentage of sequences
that must have gaps at a given position for two consecutive residues to be considered an "insert".
The second parameter "Score Treshold" controls how many sequences should be dropped based on the
score the received."""

from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import pyqtSlot
from typing import List, Optional

from ...core.Qt.QtWidgets import throttle
from ..cleanup import score_inserts
from .ScoreWithScope import ScoreByPosition, ScoreWidget
from .Ui_ScoreByLongInserts import Ui_ScoreByLongInserts

class ScoreByLongInserts(ScoreWidget):

    def __init__(self):
        super().__init__()
        self.__ui = Ui_ScoreByLongInserts()
        self.__ui.setupUi(self)
        self.__score : Optional[List[List[int]]] = None

        assert(__doc__)
        self.__ui.descriptionLabel.setText(__doc__.replace("\n", " "))
        self.__alignment : Optional[MultipleSeqAlignment] = None

        self.__ui.continuityTresholdSlider.valueChanged.connect(self.__update_score)

    @pyqtSlot(name="__update_score")
    @throttle(1000)
    def __update_score(self):

        if self.__alignment is None:
            return

        self.__score = score_inserts(
                self.__alignment,
                self.__continuity_treshold
            )

        self.on_score_changed.emit()

    @property
    def __continuity_treshold(self) -> float:
        return self.__ui.continuityTresholdSlider.value() / 100

    def score_alignment(self, alignment : MultipleSeqAlignment) -> ScoreByPosition:
        self.__alignment = alignment
        msa_score = score_inserts(alignment, self.__continuity_treshold)
        self.__score = msa_score
        return msa_score

    @property
    def score(self) -> Optional[ScoreByPosition]:
        return self.__score