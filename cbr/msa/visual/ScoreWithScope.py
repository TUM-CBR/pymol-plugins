import numpy as np
from numpy.typing import NDArray
from PyQt5.QtCore import pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import QWidget
from typing import Tuple

from ...core.Qt.QtCore import block_signals
from ...core.Qt.QtWidgets import slider_with_label
from ..cleanup.support import normalized, ranked
from .MsaCleanerResult import MsaCleanerBase, MsaCleanerResult, ScoreContext
from .Ui_ScoreWithScope import Ui_ScoreWithScope

from typing import Optional

ScoreByPosition = NDArray[np.float64]

class ScoreWidget(QWidget):

    on_score_changed = pyqtSignal()

    @property
    def score(self) -> Optional[ScoreByPosition]:
        raise NotImplementedError()

    def score_alignment(self, context : ScoreContext) -> ScoreByPosition:
        raise NotImplementedError()

class ScoreWithScope(MsaCleanerBase):

    def __init__(self, score_widget : ScoreWidget):
        super().__init__()

        self.__ui = Ui_ScoreWithScope()
        self.__ui.setupUi(self)

        self.__score_widget = score_widget
        self.__ui.gridLayout.replaceWidget(
            self.__ui.scoreWidget,
            score_widget
        )
        self.__score_widget.on_score_changed.connect(self.__on_score_by_position)

        self.__scores_by_positions : Optional[ScoreByPosition] = None
        self.__score : Optional[MsaCleanerResult] = None
        self.__treshold_slider = slider_with_label(
            self.__ui.tresholdSlider,
            self.__ui.tresholdLabel,
            1/100
        )
        self.__treshold_slider.slider_released.connect(self.__on_treshold_changed)

        self.__scope_start = slider_with_label(
            self.__ui.scopeStartSlider,
            self.__ui.scopeStartLabel,
        )
        self.__scope_start.slider_released.connect(self.__on_scope_changed)

        self.__scope_end = slider_with_label(
            self.__ui.scopeEndSlider,
            self.__ui.scopeEndLabel,
        )
        self.__scope_end.slider_released.connect(self.__on_scope_changed)

    def __update_score(self, score : Optional[MsaCleanerResult]):
        self.__score = score
        self.on_score_changed.emit()

    @pyqtSlot()
    def __on_score_by_position(self):
        self.__update_scores_by_position(self.__score_widget.score)

    def __update_scores_by_position(
        self,
        scores_by_positions: Optional[ScoreByPosition],
        force_reset_treshold: bool = False
    ):

        self.__scores_by_positions = scores_by_positions
        if scores_by_positions is None or len(scores_by_positions) == 0:
            self.__scores_by_positions = None
            self.__update_score(None)
            return

        msa_length = len(scores_by_positions[0])
        assert all(len(score) == msa_length for score in scores_by_positions)

        with block_signals(self.__ui.scopeStartSlider) \
            , block_signals(self.__ui.scopeEndSlider):
            max_scope = self.__scope_end.value
            self.__ui.scopeEndSlider.setMaximum(msa_length)
            self.__ui.scopeStartSlider.setMaximum(msa_length)

            if force_reset_treshold:
                self.__ui.scopeStartSlider.setValue(0)

            if max_scope > msa_length or force_reset_treshold:
                self.__ui.scopeEndSlider.setValue(msa_length)

        self.__scope_start.reset_label()
        self.__scope_end.reset_label()

        self.__update_score(
            MsaCleanerResult(
                self.__get_aggregate_scores(scores_by_positions),
                self.__treshold
            )
        )

    @property
    def __treshold(self) -> Tuple[float, float]:
        return (0, self.__treshold_slider.value)

    @pyqtSlot(name="__on_treshold_changed")
    def __on_treshold_changed(self):

        if self.__score is None:
            return

        self.__update_score(self.__score._replace(treshold = self.__treshold))

    def __get_aggregate_scores(self, scores: ScoreByPosition) -> NDArray[np.float64]:
        scope_start = self.__scope_start.value - 1
        scope_end = self.__scope_end.value

        result = ranked(np.sum(scores[:,scope_start:scope_end], axis=1))

        return normalized(result.astype(np.float64))

    @pyqtSlot(name="__on_scope_changed")
    def __on_scope_changed(self):

        scores_by_position = self.__scores_by_positions
        if self.__score is None or scores_by_position is None:
            return

        aggregate_scores = self.__get_aggregate_scores(scores_by_position)
        self.__update_score(
            self.__score._replace(scores = aggregate_scores)
        )

    def set_scope(self, start: int, end: int):
        self.__ui.scopeStartSlider.setValue(start)
        self.__ui.scopeEndSlider.setValue(end)
        self.__on_scope_changed()

    @property
    def score(self) -> Optional[MsaCleanerResult]:
        return self.__score

    def score_alignment(self, context : ScoreContext) -> MsaCleanerResult:

        with block_signals(self):
            score_by_alignment = self.__score_widget.score_alignment(context)
            self.__update_scores_by_position(score_by_alignment, force_reset_treshold=True)

            assert self.score is not None

            return self.score