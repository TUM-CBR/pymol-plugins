from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import QWidget
from typing import Tuple

from ...core.Qt.QtCore import block_signals
from ...core.Qt.QtWidgets import slider_with_label, throttle
from ..cleanup.support import normalized
from .MsaCleanerResult import MsaCleanerBase, MsaCleanerResult
from .Ui_ScoreWithScope import Ui_ScoreWithScope

from typing import List, Optional

ScoreByPosition = List[List[float]]

class ScoreWidget(QWidget):

    on_score_changed = pyqtSignal()

    @property
    def score(self) -> Optional[ScoreByPosition]:
        raise NotImplementedError()

    def score_alignment(self, alignment : MultipleSeqAlignment) -> ScoreByPosition:
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
        self.__treshold_slider.value_changed.connect(self.__on_treshold_changed)

        self.__scope_start = slider_with_label(
            self.__ui.scopeStartSlider,
            self.__ui.scopeStartLabel,
        )
        self.__scope_start.value_changed.connect(self.__on_scope_changed)

        self.__scope_end = slider_with_label(
            self.__ui.scopeEndSlider,
            self.__ui.scopeEndLabel,
        )
        self.__scope_end.value_changed.connect(self.__on_scope_changed)

    def __update_score(self, score : Optional[MsaCleanerResult]):
        self.__score = score
        self.on_score_changed.emit()

    @pyqtSlot()
    def __on_score_by_position(self):
        self.__update_scores_by_position(self.__score_widget.score)

    def __update_scores_by_position(self, scores_by_positions: Optional[ScoreByPosition]):

        self.__scores_by_positions = scores_by_positions
        if scores_by_positions is None or len(scores_by_positions) == 0:
            self.__scores_by_positions = None
            self.__update_score(None)
            return

        msa_length = len(scores_by_positions[0])
        assert all(len(score) == msa_length for score in scores_by_positions)

        self.__ui.scopeEndSlider.setMaximum(msa_length)
        self.__ui.scopeStartSlider.setMaximum(msa_length)
        self.__ui.scopeEndSlider.setValue(msa_length)

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
    @throttle(1000)
    def __on_treshold_changed(self):

        if self.__score is None:
            return

        self.__update_score(self.__score._replace(treshold = self.__treshold))

    def __get_aggregate_scores(self, scores: ScoreByPosition) -> List[float]:
        scope_start = self.__scope_start.value - 1
        scope_end = self.__scope_end.value

        return normalized([
            sum(score[scope_start:scope_end])
            for score in scores
        ])

    @pyqtSlot(name="__on_scope_changed")
    @throttle(1000)
    def __on_scope_changed(self):

        scores_by_position = self.__scores_by_positions
        if self.__score is None or scores_by_position is None:
            return

        aggregate_scores = self.__get_aggregate_scores(scores_by_position)
        self.__update_score(
            self.__score._replace(scores = aggregate_scores)
        )

    @property
    def score(self) -> Optional[MsaCleanerResult]:
        return self.__score

    def score_alignment(self, alignment : MultipleSeqAlignment) -> MsaCleanerResult:

        with block_signals(self):
            score_by_alignment = self.__score_widget.score_alignment(alignment)
            self.__update_scores_by_position(score_by_alignment)

            assert self.score is not None

            return self.score