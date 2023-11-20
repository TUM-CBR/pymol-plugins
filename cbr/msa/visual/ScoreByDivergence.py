"""This score is defined as the number of gaps introduced by a sequence divided by the number the largest
number of gaps instroduced by a single sequence. The higher the score means that a greater percentage of the
gaps are introduced by the sequence.
"""

from Bio.Align import MultipleSeqAlignment
from typing import cast, Optional

from ..cleanup import score_by_gap_divergence
from .Ui_ScoreByDivergence import Ui_ScoreByDivergence
from .ScoreWithScope import ScoreByPosition, ScoreWidget

class ScoreByDivergence(ScoreWidget):

    def __init__(self):
        super().__init__()
        self.__ui = Ui_ScoreByDivergence()
        self.__ui.setupUi(self)
        assert(__doc__)
        self.__ui.descriptionLabel.setText(__doc__.replace("\n", " "))
        self.__score : Optional[ScoreByPosition] = None

    @property
    def score(self) -> Optional[ScoreByPosition]:
        return self.__score

    def score_alignment(self, alignment : MultipleSeqAlignment) -> ScoreByPosition:
        msa_score = score_by_gap_divergence(alignment)
        score = cast(ScoreByPosition, msa_score)
        self.__score = score
        return score
