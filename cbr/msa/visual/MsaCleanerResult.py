import numpy as np
from numpy.typing import NDArray
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QWidget
from typing import NamedTuple, Optional, Tuple

from ..cleanup.support import ScoreContext

class MsaCleanerResult(NamedTuple):
    scores : NDArray[np.float64]
    treshold : Tuple[float, float]

class MsaCleanerBase(QWidget):

    on_score_changed = pyqtSignal()

    @property
    def score(self) -> Optional[MsaCleanerResult]:
        raise NotImplemented()

    def score_alignment(self, context : ScoreContext) -> MsaCleanerResult:
        raise NotImplemented()

    def set_scope(self, start: int, end: int):
        pass