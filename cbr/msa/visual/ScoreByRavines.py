"""This scoring function is designed to eliminate sequences produce ravines in a multiple
sequence alignment. A ravine is a long section where most sequences contain gaps. This is
sometimes created by a few sequences having an insert in that section but they may also
be created by multiple sequences having a single residue in several of the positions. The
'Continuity Treshold' parameter controls the percentage of sequences that must contain a
gap at a given position for that position to be considered a gap. For example, a value of
'0.9' means that 90% of the sequences must contain a gap at a position for that position
to be marked as a gap."""

from typing import Optional

from PyQt5.QtCore import pyqtSlot
import numpy as np
from numpy.typing import NDArray

from ...control import viter
from ...core.Qt.QtCore import block_signals
from ...core.Qt.QtWidgets import slider_with_label, throttle
from ..cleanup import score_residues_in_ravines
from ..cleanup.support import ScoreContext
from .ScoreWithScope import ScoreByPosition, ScoreWidget
from .Ui_ScoreByRavines import Ui_ScoreByRavines

class ScoreByRavines(ScoreWidget):
    """This widget is responsible for running the :func:`score_residues_in_ravines`
    scoring on a Multiple Sequence Alignment. It provides a slider widget that allows
    the user to configure the 'continuity_treshold' parameter of the function."""

    def __init__(self):
        super().__init__()
        self.__ui = Ui_ScoreByRavines()
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
        self.__continuity_treshold_slider.value_changed.connect(self.__on_continuity_treshold_changed)

    def __set_score(self, score: NDArray[np.int64]) -> ScoreByPosition:
        """Updates the `ScoreByRavines.__score` from a scoring with int rather
        than float as that is what the :func:`score_residues_in_ravines` returns.
        This is mostly there to maye the typechecker happy."""
        self.__score = result = score.astype(np.float64)
        self.on_score_changed.emit()
        return result

    @property
    def __continuity_treshold(self) -> float:
        return self.__continuity_treshold_slider.value

    @pyqtSlot(name="__on_continuity_treshold_changed")
    @throttle(1000)
    def __on_continuity_treshold_changed(self):
        """This function runs the :func:`score_residues_in_ravines` for the current
        alignment to get an updated score. It will not do anything if there is no
        alignment as that means that the user is playing with the controls w/o having
        any alignment loaded."""

        for alignment in viter(self.__alignment):
            self.__set_score(
                score_residues_in_ravines(
                    alignment,
                    self.__continuity_treshold
                )
            )

    def score_alignment(self, context : ScoreContext) -> ScoreByPosition:
        """Update the multiple sequence alignment being scored and run the initial
        scoring. This function will not trigger any signals as the result is intended
        to be the first score of the new alignment."""

        with block_signals(super()):
            self.__alignment = context
            msa_score = score_residues_in_ravines(context, self.__continuity_treshold)
            return self.__set_score(msa_score)

    @property
    def score(self) -> Optional[ScoreByPosition]:
        return self.__score

    