"""This scoring function is designed to eliminate sequences that contain long inserts in the MSA.
It has two parameters. The first one "Continuity Treshold" determines the percentage of sequences
that must have gaps at a given position for two consecutive residues to be considered an "insert".
The second parameter "Score Treshold" controls how many sequences should be dropped based on the
score the received."""

from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import pyqtSlot
from typing import cast, List, Optional

from ...core.Qt.QtWidgets import throttle
from ..cleanup import score_inserts
from .ScoreWithScope import ScoreByPosition, ScoreWidget
from .Ui_ScoreByLongInserts import Ui_ScoreByLongInserts

class ScoreByLongInserts(ScoreWidget):
    """This widget is responsible for running the :func:`score_inserts` scoring
    function on a Multiple Sequence Alignment. It provides a slider that allows
    the user to configure the 'continuity_treshold' parameter of the function"""

    def __init__(self):
        super().__init__()
        self.__ui = Ui_ScoreByLongInserts()
        self.__ui.setupUi(self)
        self.__score : Optional[ScoreByPosition] = None

        assert(__doc__)
        self.__ui.descriptionLabel.setText(__doc__.replace("\n", " "))
        self.__alignment : Optional[MultipleSeqAlignment] = None

        self.__ui.continuityTresholdSlider.valueChanged.connect(self.__update_score)

    def __set_score(self, score: List[List[int]]) -> ScoreByPosition:
        value = cast(ScoreByPosition, score)
        self.__score = value
        return value

    @pyqtSlot(name="__update_score")
    @throttle(1000)
    def __update_score(self):

        if self.__alignment is None:
            return

        self.__set_score(
            score_inserts(
                self.__alignment,
                self.__continuity_treshold
            )
        )

        self.on_score_changed.emit()

    @property
    def __continuity_treshold(self) -> float:
        return self.__ui.continuityTresholdSlider.value() / 100

    def score_alignment(self, alignment : MultipleSeqAlignment) -> ScoreByPosition:

        self.__alignment = alignment
        msa_score = score_inserts(alignment, self.__continuity_treshold)
        return self.__set_score(msa_score)

    @property
    def score(self) -> Optional[ScoreByPosition]:
        return self.__score