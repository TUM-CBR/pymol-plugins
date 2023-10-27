from Bio.Align import MultipleSeqAlignment
from typing import Optional, Tuple

from ..cleanup import score_by_gap_divergence
from .Ui_ScoreByDivergence import Ui_ScoreByDivergence
from .MsaCleanerResult import MsaCleanerBase, MsaCleanerResult

class ScoreByDivergence(MsaCleanerBase):

    def __init__(self):
        super().__init__()
        self.__ui = Ui_ScoreByDivergence()
        self.__ui.setupUi(self)
        self.__score : Optional[MsaCleanerResult] = None

    @property
    def score(self) -> Optional[MsaCleanerResult]:
        return self.__score

    @property
    def __treshold(self) -> Tuple[float, float]:
        return (0, self.__ui.divergenceSlider.value() / 100)

    def score_alignment(self, alignment : MultipleSeqAlignment) -> MsaCleanerResult:
        msa_score = score_by_gap_divergence(alignment)
        result = MsaCleanerResult(
            scores=msa_score,
            treshold=self.__treshold
        )

        self.__score = result
        return result
