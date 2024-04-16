from operator import is_
from statistics import mean
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from numpy.typing import NDArray
import pymol
from PyQt5.QtCore import pyqtSlot, QItemSelection, QModelIndex, QObject, Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QAbstractItemView, QDialog, QWidget
from typing import Any, cast, Dict, List, NamedTuple, Optional, Sequence

from .AbstractStructureTableModel import BLACK, AbstractStructureTableModel, StructureMapping, StructureMappings
from .MultiStructureSelector import MultiStructureSelector
from .Ui_StructureCompare import Ui_StructureCompare

from ...core.Context import Context
from ...core.pymol.algorithms import cx_coords
from ...core.pymol.structure import StructureSelection
from ...clustal.msa import is_blank
from ...clustal.Clustal import Clustal, get_clustal_from_context

from ..display.sequence import RESIDUE_COLORS

def score_random(msa: MultipleSeqAlignment) -> List[Optional[float]]:
    from random import random

    return [
        2*random() - 1
        for _ in range(0, msa.get_alignment_length())
    ]

def score_fragment_identity(msa: MultipleSeqAlignment) -> List[Optional[float]]:

    i = 0
    end = msa.get_alignment_length()
    result: List[Optional[float]] = list(None for _ in range(end))

    assert len(msa) == 2, "Score fragment identity requries an alignment of two sequences"

    seq_1: Sequence[str] = cast(Sequence[str], msa[0]._seq)
    seq_2: Sequence[str] = cast(Sequence[str], msa[1]._seq)

    while i < end:

        r1 = seq_1[i].upper()
        r2 = seq_2[i].upper()
        start = i
        matches = 0
        while i < end and not is_blank(r1) and not is_blank(r2):
            matches += 1 if r1 == r2 else 0
            i += 1

            if i < end:
                r1 = seq_1[i].upper()
                r2 = seq_2[i].upper()

        span = i - start
        identity = matches / span if span > 0 else 0
        for j in range(start, i):
            result[j] = identity

        if span == 0:
            i += 1

    return result

def score_residue_distance(
    msa: MultipleSeqAlignment,
    resv_mapping_1: List[int],
    resv_mapping_2: List[int],
    structure_1: StructureSelection,
    structure_2: StructureSelection
) -> List[Optional[float]]:
    assert len(msa) == 2, "Score fragment identity requries an alignment of two sequences"

    end = msa.get_alignment_length()
    result : List[Optional[float]] = [None for _ in range(end)]
    seq_mapping_1 = msa.alignment.indices[0]
    seq_mapping_2 = msa.alignment.indices[1]
    coords_1 = cx_coords(structure_1).resv_to_position(0)
    coords_2 = cx_coords(structure_2).resv_to_position(0)

    for i in range(end):

        pos_1 = seq_mapping_1[i]
        pos_2 = seq_mapping_2[i]

        if pos_1 < 0 or pos_2 < 0:
            continue

        resv_1 = resv_mapping_1[pos_1]
        resv_2 = resv_mapping_2[pos_2]

        pos_1 = coords_1.get(resv_1)
        pos_2 = coords_2.get(resv_2)

        # Residues may appear in the sequence but not be visisble
        if pos_1 is None or pos_2 is None:
            continue

        result[i] = float(np.linalg.norm(pos_2 - pos_1))

    return result

COLOR_RANGE = 255

class ScoreEntry(NamedTuple):
    absolute: float
    relative: float

    def to_qcolor(self) -> QColor:

        score = self.relative

        return QColor(
            max(0, min(255, int(-166*(score - 0.5)))),
            max(0, min(255, int(166*(score + 0.5)))),
            0
        )

MsaScore = List[Optional[ScoreEntry]]

def to_score_entries(scores: List[Optional[float]], invert_ranking: bool = False) -> MsaScore:
    only_values = [s for s in scores if s is not None]
    max_score = max(only_values)
    min_score = min(only_values)
    avg_score = mean(only_values)
    pos_range = max_score - avg_score
    neg_range = avg_score - min_score
    range_sgn = -1 if invert_ranking else 1

    def get_entry(score: Optional[float]) -> Optional[ScoreEntry]:

        if score is None:
            return None
        elif score > avg_score:
            rel_score = abs((score - avg_score) / pos_range)
        elif score < avg_score:
            rel_score = -abs((avg_score - score) / neg_range)
        else:
            rel_score = 0

        return ScoreEntry(score, range_sgn*rel_score)

    return [
        get_entry(score)
        for score in scores
    ]

class DataEntry(NamedTuple):
    sequence: Optional[List[str]] = None
    scores: Optional[MsaScore] = None

K_SCORE_IDENTITY = "Identity"
K_SCORE_DISTANCE = "Distance"
K_SCORE_NONE = "None"
K_SCORES = [K_SCORE_IDENTITY, K_SCORE_DISTANCE, K_SCORE_NONE]

class StructureMsaEntry(NamedTuple):
    structures: List[StructureSelection]
    resv_mappings: List[List[int]]
    msa: MultipleSeqAlignment
    structure_mappings: StructureMappings
    scores: Dict[str, MsaScore]
    data: Dict[str, DataEntry]

    @classmethod
    def to_structure_mappings(cls, msa: MultipleSeqAlignment, i: int, resv_mappings: List[int]) -> List[int]:
        r_indexes: NDArray[np.int64] = msa.alignment.indices[i]
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
            K_SCORE_IDENTITY: to_score_entries(score_fragment_identity(msa)),
            K_SCORE_DISTANCE: to_score_entries(
                score_residue_distance(
                    msa,
                    resv_mappings[0],
                    resv_mappings[1],
                    structures[0],
                    structures[1]
                ),
                invert_ranking=True
            )
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
        selected_score: str,
        parent: Optional[QObject] = None,
    ) -> None:
        super().__init__(parent)
        self.__clustal = clustal
        self.__entry: Optional[StructureMsaEntry] = None
        self.__selected_score: str = selected_score
        self.__rows: List[str] = []

    def set_structures(
        self,
        structure_a: StructureSelection,
        structure_b: StructureSelection
    ):
        pymol.cmd.align(
            structure_a.selection,
            structure_b.selection
        )

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

        self.__rows = \
            [msa[0].id] + \
            list(self.__entry.scores.keys()) + \
            [msa[1].id]

        self.modelReset.emit()
        self.structure_mappings_signal.emit()

    def set_selected_score(self, score: str):
        self.__selected_score = score
        self.structure_mappings_signal.emit()

    def __get_structure_color__(self, row: Optional[int], col: Optional[int]) -> Optional[QColor]:
        entry = self.__entry
        selected = self.__selected_score
        if col is None \
            or entry is None \
            or selected is None \
            or selected not in entry.scores:
            return None
        
        score = entry.scores[selected][col]

        if score is None:
            return None
        else:
            return score.to_qcolor()

    def __get_structure_mappings__(self) -> Optional[StructureMappings]:
        entry = self.__entry

        if entry is None or self.__selected_score == K_SCORE_NONE:
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
        
    def __get_record(self, i_row: int) -> Optional[DataEntry]:

        entry = self.__entry

        if entry is None:
            return None
        row = self.__rows[i_row]
        return entry.data[row]
    
    def __get_sequence_column(self, col: int, seq: List[str]) -> str:
        return seq[col]
    
    def __get_score(self, col: int, scores: MsaScore):
        return scores[col]

    def __display_role_data(self, index: QModelIndex) -> Optional[str]:

        entry = self.__get_record(index.row())

        if entry is None:
            return None
        elif entry.sequence:
            return self.__get_sequence_column(index.column(), entry.sequence)
        elif entry.scores:
            score = self.__get_score(index.column(), entry.scores)
            return str(round(score.absolute, 2)) if score is not None else None
        else:
            return None
        
    def __background_display_role(self, index: QModelIndex) -> Optional[QColor]:

        entry = self.__get_record(index.row())

        if entry is None:
            return None
        elif entry.sequence:
            value = self.__get_sequence_column(index.column(), entry.sequence).upper()
            return RESIDUE_COLORS.get(value)
        elif entry.scores:
            score = self.__get_score(index.column(), entry.scores)
            return BLACK if score is None else score.to_qcolor()
        else:
            return None

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if role == Qt.ItemDataRole.DisplayRole:
            return self.__display_role_data(index)
        elif role == Qt.ItemDataRole.BackgroundColorRole:
            return self.__background_display_role(index)
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

        clustal: Clustal = get_clustal_from_context(context)
        self.__model = StructureCompareModel(clustal, K_SCORES[0])
        self.__select_structures = MultiStructureSelector(2,2)
        self.__ui.selectStructuresButton.clicked.connect(self.__on_select_structures_clicked)
        table = self.__ui.structuresTable

        table.setModel(self.__model)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectColumns)

        selection_model = table.selectionModel()
        assert  selection_model is not None, "Selection model unexpectedly None"

        selection_model.selectionChanged.connect(self.__on_selection_changed)

        structures_combo = self.__ui.structureColoringCombo
        structures_combo.addItems(K_SCORES)
        structures_combo.setCurrentIndex(0)
        structures_combo.currentIndexChanged.connect(self.__on_score_selected)

    @pyqtSlot(int)
    def __on_score_selected(self, index: int):
        if index < 0:
            return
        
        self.__model.set_selected_score(K_SCORES[index])

    @pyqtSlot(QItemSelection, QItemSelection)
    def __on_selection_changed(self, selected: QItemSelection, deselected: QItemSelection):
        self.__model.update_selection(
            self.__ui.structuresTable.selectedIndexes()
        )

    @pyqtSlot()
    def __on_select_structures_clicked(self):

        result = self.__select_structures.exec()

        if result == QDialog.DialogCode.Accepted:
            structures = self.__select_structures.selected_values()
            self.__model.set_structures(
                structures[0],
                structures[1]
            )
            self.__ui.structuresTable.resizeColumnsToContents()



