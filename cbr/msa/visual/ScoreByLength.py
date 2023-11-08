from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import pyqtSlot
from typing import List, Optional, Tuple

from ...core.Qt.QtWidgets import throttle
from ..cleanup import score_length
from .MsaCleanerResult import MsaCleanerBase, MsaCleanerResult
from .Ui_ScoreByLength import Ui_ScoreByLength

class ScoreByLength(MsaCleanerBase):

    def __init__(self):
        super().__init__()
        self.__ui = Ui_ScoreByLength()
        self.__ui.setupUi(self)
        self.__ui.minSlider.valueChanged.connect(self.__on_slider_value_changed)
        self.__ui.maxSlider.valueChanged.connect(self.__on_slider_value_changed)
        self.__score : Optional[MsaCleanerResult] = None
        self.__score_widgets = [
            (self.__ui.minLengthLabel, self.__ui.minSlider),
            (self.__ui.maxLenghtLabel, self.__ui.maxSlider)
        ]

    @pyqtSlot()
    def __on_slider_value_changed(self):

        if self.__score is None:
            return

        for (label, slider) in self.__score_widgets:
            label.setText(str(slider.value()))

        self.__score = self.__score._replace(treshold = self.__treshold)
        self.__on_min_value_changed()

    @property
    def __treshold(self) -> Tuple[int, int]:
        return (
            self.__ui.minSlider.value(),
            self.__ui.maxSlider.value()
        )

    @throttle(1000)
    def __on_min_value_changed(self):
        self.on_score_changed.emit()

    @property
    def score(self) -> Optional[MsaCleanerResult]:
        return self.__score

    def __set_range(self, scores: List[float]) -> None:
        min_range = int(min(scores))
        max_range = int(max(scores))

        for (_, slider) in self.__score_widgets:
            slider.setMinimum(min_range)
            slider.setMaximum(max_range)

        self.__ui.maxSlider.setValue(max_range)
        
        for (label, slider) in self.__score_widgets:
            label.setText(str(slider.value()))

    def score_alignment(self, alignment : MultipleSeqAlignment) -> MsaCleanerResult:
        msa_score = score_length(alignment)

        self.__set_range(msa_score)

        result = MsaCleanerResult(
            scores=msa_score,
            treshold=self.__treshold
        )

        self.__score = result
        return result