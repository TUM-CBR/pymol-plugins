"""This scoring function is designed to eliminate sequences that contain long inserts in the MSA.
It has two parameters. The first one "Continuity Treshold" determines the percentage of sequences
that must have gaps at a given position for two consecutive residues to be considered an "insert".
"""
import numpy as np
from numpy.typing import NDArray
from PyQt5.QtCore import pyqtSlot
from typing import cast, Optional

from ...core.Qt.QtWidgets import slider_with_label, throttle
from ..cleanup import score_inserts
from .ScoreWithScope import ScoreByPosition, ScoreContext, ScoreWidget
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
        self.__alignment : Optional[ScoreContext] = None
        self.__continuity_treshold_slider = slider_with_label(
            self.__ui.continuityTresholdSlider,
            self.__ui.continuityTresholdLabel,
            factor=1/100
        )

        self.__continuity_treshold_slider.value_changed.connect(self.__update_score)

    def __set_score(self, score: NDArray[np.int64]) -> ScoreByPosition:
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
        return self.__continuity_treshold_slider.value

    def score_alignment(self, context : ScoreContext) -> ScoreByPosition:

        self.__alignment = context
        msa_score = score_inserts(context, self.__continuity_treshold)
        return self.__set_score(msa_score)

    @property
    def score(self) -> Optional[ScoreByPosition]:
        return self.__score