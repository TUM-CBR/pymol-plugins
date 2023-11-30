from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import QRect, pyqtSlot
from typing import List, Optional, Tuple

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPaintEvent, QPixmap
from PyQt5.QtWidgets import QWidget

from ...core.Qt.QtGui import paint
from ...core.Qt.QtWidgets import throttle
from ..cleanup import score_length
from .MsaCleanerResult import MsaCleanerBase, MsaCleanerResult
from .Ui_ScoreByLength import Ui_ScoreByLength

class FreqCanvas(QWidget):

    def __init__(self):
        super().__init__()
        self.__pixmap = None
        self.setMinimumHeight(100)

    def update_pixmap(self, pixmap : Optional[QPixmap]):
        self.__pixmap = pixmap
        self.repaint()

    def paintEvent(self, a0: QPaintEvent) -> None:
        super().paintEvent(a0)

        if self.__pixmap is not None:
            with paint(self) as ctx:
                ctx.painter.drawPixmap(
                    QRect(0,0, self.width(), self.height()),
                    self.__pixmap
                )

class ScoreByLength(MsaCleanerBase):

    def __init__(self):
        super().__init__()
        self.__ui = Ui_ScoreByLength()
        self.__ui.setupUi(self)
        self.__ui.minSlider.valueChanged.connect(self.__on_slider_value_changed)
        self.__ui.maxSlider.valueChanged.connect(self.__on_slider_value_changed)
        self.__score : Optional[MsaCleanerResult] = None
        self.__freq_canvas = FreqCanvas()
        self.layout().replaceWidget(
            self.__ui.frequencyWidget,
            self.__freq_canvas
        )
        self.__freq_canvas.show()
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

    def __render_length_frequency(self, scores: List[float]):
        min_range = int(min(scores))
        max_range = int(max(scores))
        score_range = max_range - min_range

        # On average, we want 20 sequences
        # per bucket
        buckets = 1 + round(len(scores) / 20)

        if buckets < 10:
            buckets = 10

        bucket_size = score_range / buckets
        freq = [0 for _ in range(0,buckets)]

        for score in scores:
            bucket = int(min((score - min_range) / bucket_size, buckets - 1))
            freq[bucket] += 1

        pixmap_h = max(freq)

        pixmap = QPixmap(buckets, pixmap_h)
        with paint(pixmap) as context:
            painter = context.painter

            painter.setPen(Qt.white)
            painter.setBrush(Qt.white)
            painter.drawRect(0,0, buckets, pixmap_h)

            painter.setPen(Qt.darkBlue)
            painter.setBrush(Qt.darkBlue)

            for i,count in enumerate(freq):
                painter.drawLine(i, pixmap_h - count, i, pixmap_h)

        return pixmap

    def __set_range(self, scores: List[float]) -> None:
        min_range = int(min(scores))
        max_range = int(max(scores))

        self.__freq_canvas.update_pixmap(self.__render_length_frequency(scores))

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

