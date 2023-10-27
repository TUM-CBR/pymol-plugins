from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QWidget
from typing import List, NamedTuple, Tuple

class MsaCleanerResult(NamedTuple):
    scores : List[float]
    treshold : Tuple[float, float]

class MsaCleanerResultImplementation(MsaCleanerResult):

    def __init__(
        self,
        scores : List[float]
    ):
        self.__scores = scores

    @property
    def scores(self) -> List[float]:
        return self.__scores

class MsaCleanerBase(QWidget):

    on_score_changed = pyqtSignal()

    def score_alignment(self, alignment : MultipleSeqAlignment) -> MsaCleanerResult:
        raise NotImplemented()