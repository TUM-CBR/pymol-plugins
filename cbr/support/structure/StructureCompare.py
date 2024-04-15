from statistics import mean
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from PyQt5.QtGui import QColor
import numpy as np
from numpy.typing import NDArray
from PyQt5.QtCore import QModelIndex, QObject, Qt
from PyQt5.QtWidgets import QWidget
from typing import Any, Dict, List, NamedTuple, Optional

from .AbstractStructureTableModel import AbstractStructureTableModel, StructureMapping, StructureMappings
from .Ui_StructureCompare import Ui_StructureCompare

from ...core.Context import Context
from ...core.pymol.structure import StructureSelection
from ...clustal.Clustal import Clustal, get_clustal_from_context

def score_random(msa: MultipleSeqAlignment) -> List[Optional[float]]:
    from random import random

    return [
        2*random() - 1
        for _ in range(0, msa.get_alignment_length())
    ]

class ScoreEntry(NamedTuple):
    absolute: float
    relative: float

    def to_qcolor(self) -> QColor:

        score = self.relative

        if score == 0:
            return QColor(255,255,255)
        elif score > 0:
            return QColor(0, min(255, int(score * 255)), 0)
        else:
            score = abs(score)
            return QColor(min(255, int(score * 255)), 0, 0)


MsaScore = List[Optional[ScoreEntry]]

def to_score_entries(scores: List[Optional[float]]) -> MsaScore:
    only_values = [s for s in scores if s is not None]
    max_score = max(only_values)
    min_score = min(only_values)
    avg_score = mean(only_values)
    pos_range = max_score - avg_score
    neg_range = avg_score = min_score

    def get_entry(score: Optional[float]) -> Optional[ScoreEntry]:

        if score is None:
            return None
        elif score > avg_score:
            rel_score = abs((score - avg_score) / pos_range)
        elif score < avg_score:
            rel_score = - abs((avg_score - score) / neg_range)
        else:
            rel_score = 0

        return ScoreEntry(score, rel_score)

    return [
        get_entry(score)
        for score in scores
    ]

class DataEntry(NamedTuple):
    sequence: Optional[List[str]] = None
    scores: Optional[MsaScore] = None

class StructureMsaEntry(NamedTuple):
    structures: List[StructureSelection]
    resv_mappings: List[List[int]]
    msa: MultipleSeqAlignment
    structure_mappings: StructureMappings
    scores: Dict[str, MsaScore]
    data: Dict[str, DataEntry]

    @classmethod
    def to_structure_mappings(cls, msa: MultipleSeqAlignment, i: int, resv_mappings: List[int]) -> List[int]:
        r_indexes: NDArray[np.int64] = msa.alignment.inverse_indices[i]
        last: int = resv_mappings[0]

        def get_resv(i: int) -> int:
            nonlocal r_indexes
            nonlocal last
            index = r_indexes[i]

            if i >= 0:
                last = resv_mappings[index]

            return last

        return [
            get_resv(i)
            for i in range(0, msa.get_alignment_length())
        ]

    @classmethod
    def to_structures_mappings(
        cls,
        msa: MultipleSeqAlignment,
        structures: List[StructureSelection],
        resv_mappings: List[List[int]]
    ) -> StructureMappings:

        return StructureMappings(
            mappings=[
                StructureMapping(
                    struct,
                    by_column=cls.to_structure_mappings(msa, i, resv_mappings[i])
                )
                for i, struct in enumerate(structures)
            ]
        )
    
    @classmethod
    def create(
        cls,
        structures: List[StructureSelection],
        resv_mappings: List[List[int]],
        msa: MultipleSeqAlignment
    ):
        
        scores = {
            "random 1": to_score_entries(score_random(msa)),
            "random 2": to_score_entries(score_random(msa))
        }

        data: Dict[str, DataEntry] = {}

        for name, score in scores.items():
            data[name] = DataEntry(scores = score)

        for seq in msa:
            data[seq.id] = DataEntry(sequence = list(seq._seq))
        
        return StructureMsaEntry(
            structures,
            resv_mappings,
            msa,
            cls.to_structures_mappings(msa, structures, resv_mappings),
            scores,
            data
        )

class StructureCompareModel(AbstractStructureTableModel):

    def __init__(
        self,
        clustal: Clustal,
        parent: Optional[QObject] = None
    ) -> None:
        super().__init__(parent)
        self.__clustal = clustal
        self.__entry: Optional[StructureMsaEntry] = None
        self.__selected_score: Optional[str] = None
        self.__rows: List[str] = []

    def set_structures(
        self,
        structure_a: StructureSelection,
        structure_b: StructureSelection
    ):
        seq_a = structure_a.get_sequence()
        resv_a = list(seq_a.keys())
        resv_a.sort()
        
        seq_b = structure_b.get_sequence()
        resv_b = list(seq_b.keys())
        resv_b.sort()

        msa = self.__clustal.run_msa_seqs([
            SeqRecord(Seq("".join(seq_a[i] for i in resv_a)), id=structure_a.show()),
            SeqRecord(Seq("".join(seq_b[i] for i in resv_b)), id=structure_b.show())
        ])

        self.__entry = StructureMsaEntry.create(
            [structure_a, structure_b],
            [resv_a, resv_b],
            msa
        )
        self.__selected_score = next(self.__entry.scores.keys().__iter__())

        self.__rows = \
            [msa[0].id] + \
            list(self.__entry.scores.keys()) + \
            [msa[1].id]

        self.modelReset.emit()
        self.structure_mappings_signal.emit()

    def __get_structure_color__(self, row: Optional[int], col: Optional[int]) -> Optional[QColor]:
        entry = self.__entry
        selected = self.__selected_score
        if col is None or entry is None or selected is None:
            return None
        
        score = entry.scores[selected][col]

        if score is None:
            return None
        else:
            return score.to_qcolor()

    def __get_structure_mappings__(self) -> Optional[StructureMappings]:
        entry = self.__entry

        if entry is None:
            return None
        else:
            return entry.structure_mappings
        
    def rowCount(self, parent: QModelIndex = ...) -> int:
        entry = self.__entry

        if entry is None:
            return 0

        return len(self.__rows)
    
    def columnCount(self, parent: QModelIndex = ...) -> int:
        entry = self.__entry

        if entry is None:
            return 0
        
        return entry.msa.get_alignment_length()
    
    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        entry = self.__entry
        if orientation == Qt.Orientation.Vertical and role == Qt.ItemDataRole.DisplayRole and entry is not None:
            return self.__rows[section]
        else:
            return None

class StructureCompare(QWidget):

    def __init__(
        self,
        context: Context,
        parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(parent)

        self.__ui = Ui_StructureCompare()
        self.__ui.setupUi(self)

        self.__clustal: Clustal = get_clustal_from_context(context)



