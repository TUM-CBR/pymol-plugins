from Bio.Align import MultipleSeqAlignment, SeqRecord
from enum import Enum
import numpy as np
from numpy.typing import NDArray
from PyQt5.QtCore import pyqtSlot, QAbstractTableModel, QModelIndex, Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QFileDialog, QWidget
from typing import Any, Callable, Dict, Iterable, Optional, cast, List, NamedTuple, Tuple

from ...core.Context import (Context)
from ...core.Qt.QtCore import DictionaryModel
from ...core.Qt.QtWidgets import show_info, throttle, with_error_handler
from ..cleanup import ScoreContext
from ...support import msa
from ...support.msa.io import save_msa
from ...support.msa.visual.MsaViewer import MsaViewer

from .Ui_MsaCleaner import Ui_MsaCleaner
from .MsaCleanerResult import MsaCleanerBase, MsaCleanerResult
from .ScoreByDivergence import ScoreByDivergence
from .ScoreByLongInserts import ScoreByLongInserts
from .ScoreByLength import ScoreByLength
from .ScoreByRavines import ScoreByRavines
from .ScoreWithScope import ScoreWithScope

class ScoreEntry(NamedTuple):
    name : str
    scores : NDArray[np.float64]
    treshold : Tuple[float, float]

    def is_included(self, i :int) -> bool:
        low, high = self.treshold
        score = self.scores[i]
        return score >= low and score <= high

    @staticmethod
    def from_result(name : str, result : MsaCleanerResult) -> 'ScoreEntry':
        return ScoreEntry(
            name,
            result.scores,
            result.treshold
        )

    def score(self, i : int) -> float:
        return self.scores[i]

class ScoreOverride(Enum):
    Green = 0
    Red = 1
    NoOverride = 2

class ScoreMeta(NamedTuple):
    override : ScoreOverride

    def check_state(self, override : ScoreOverride) -> Qt.CheckState:
        if self.override == override:
            return Qt.CheckState.Checked
        else:
            return Qt.CheckState.Unchecked

    def toggle(self, base_status : ScoreOverride) -> 'ScoreMeta':

        if self.override != base_status:
            return self._replace(override = base_status)
        else:
            return self._replace(override = ScoreOverride.NoOverride)

class SequenceScoresModel(QAbstractTableModel):

    fixed_headers = ["Greenlist", "Redlist", "Sequence Id"]
    GREENLIST_COLUMN = fixed_headers.index("Greenlist")
    REDLIST_COLUMN = fixed_headers.index("Redlist")
    OVERRIDE_COLUMNS = [GREENLIST_COLUMN, REDLIST_COLUMN]

    def __init__(
        self,
        alignment : MultipleSeqAlignment,
        scores : List[ScoreEntry],
        meta_override: Optional[Callable[[int, ScoreMeta], ScoreMeta]] = None
    ):
        super().__init__()

        self.__alignment = alignment
        self.__scores = scores

        meta_override = meta_override if meta_override is not None else lambda _,score: score

        self.__score_meta = [
            meta_override(i, ScoreMeta(override=ScoreOverride.NoOverride))
            for i in range(0, len(alignment))
        ]

    @property
    def scores(self) -> List[ScoreEntry]:
        return self.__scores

    @property
    def meta(self) -> List[ScoreMeta]:
        return self.__score_meta

    def get_index_by_name(self) -> Dict[str, int]:
        return dict(
            (seq.id, i)
            for (i, seq) in enumerate(self.__alignment)
        )

    @property
    def alignment(self) -> MultipleSeqAlignment:
        return self.__alignment

    @property
    def headers(self):
        return self.fixed_headers + [score.name for score in self.__scores]

    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__alignment)

    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.headers)

    def set_always_keep(self, index : QModelIndex, keep : bool):
        row = index.row()
        self.__score_meta[row] = self.__score_meta[row]._replace(keep_always = keep)
        self.dataChanged.emit(index, index)

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            headers = self.headers
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def __get_seq(self, index : int) -> SeqRecord:
        return cast(SeqRecord, self.__alignment[index])

    @staticmethod
    def fromat_score(value : float) -> str:
        return str(round(value, 2))

    def toggle_override(self, index : QModelIndex, override: ScoreOverride) -> None:
        """Toggle wether to alwasy include the current
        sequence in the MSA regardless of the analysis results.
        Only toggle if the index corresponds to the correct column."""

        assert index.column() in self.OVERRIDE_COLUMNS

        self.__score_meta[index.row()] = self.__score_meta[index.row()].toggle(override)
        self.dataChanged.emit(
            index,
            index.siblingAtColumn(self.columnCount() - 1),
            [Qt.ItemDataRole.CheckStateRole, Qt.ItemDataRole.BackgroundColorRole]
        )

    def __is_included(self, index : QModelIndex) -> bool:
        return self.is_included(index.row())

    def is_included(self, row_ix : int) -> bool:
        override = self.__score_meta[row_ix].override

        return override != ScoreOverride.Red \
            and (
                override == ScoreOverride.Green \
                or all(score.is_included(row_ix) for score in self.__scores)
            )

    def update_scores(self, scores : List[ScoreEntry]) -> None:
        self.__scores = scores
        
        self.dataChanged.emit(
            self.index(0, 0),
            self.index(self.rowCount() - 1, self.columnCount() - 1)
        )

    def flags(self, index: QModelIndex) -> Qt.ItemFlags:

        if index.column() in self.OVERRIDE_COLUMNS:
            return Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsUserCheckable
        else:
            return super().flags(index)

    def setData(self, index: QModelIndex, value: Any, role: int = Qt.ItemDataRole.DisplayRole) -> bool:

        if index.column() == self.GREENLIST_COLUMN and role == Qt.ItemDataRole.CheckStateRole:
            self.toggle_override(index, ScoreOverride.Green)
            return True

        if index.column() == self.REDLIST_COLUMN and role == Qt.ItemDataRole.CheckStateRole:
            self.toggle_override(index, ScoreOverride.Red)
            return True

        return super().setData(index, value, role)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole):

        if not index.isValid():
            return None

        row_ix = index.row()
        col_ix = index.column()
        sequence = self.__get_seq(row_ix)
        if role == Qt.ItemDataRole.DisplayRole:
            cols : List[Optional[str]] = [
                "",
                "",
                sequence.id
            ] + [self.fromat_score(score.score(row_ix)) for score in self.__scores]
            return cols[col_ix]
        elif role == Qt.ItemDataRole.BackgroundColorRole:
            if self.__is_included(index):
                return None
            else:
                return QColor(255,0,0,100)
        elif role == Qt.ItemDataRole.CheckStateRole and col_ix == self.GREENLIST_COLUMN:
            return self.__score_meta[row_ix].check_state(ScoreOverride.Green)
        elif role == Qt.ItemDataRole.CheckStateRole and col_ix == self.REDLIST_COLUMN:
            return self.__score_meta[row_ix].check_state(ScoreOverride.Red)
        else:
            return None


class MsaCleaner(QWidget):

    def __init__(self, ctx : Context):
        super().__init__()
        self.__ui = Ui_MsaCleaner()
        self.__ui.setupUi(self)
        self.__scores_model : Optional[SequenceScoresModel] = None

        self.__cleaners : List[Tuple[str, MsaCleanerBase]] = [
            ("Gap Divergence", ScoreWithScope(ScoreByDivergence())),
            ("Eliminate Long Inserts", ScoreWithScope(ScoreByLongInserts())),
            ("Eliminate Large Gaps", ScoreWithScope(ScoreByRavines())),
            ("Sequence Length", ScoreByLength())
        ]

        for name, cleaner in self.__cleaners:
            self.__ui.cleanersWidget.addTab(cleaner, name)
            cleaner.on_score_changed.connect(self.__on_score_changed)

        self.__msa_selector = msa.msa_selector(self, self.__ui.loadMsaButton, self.__ui.selectedFileLabel)
        self.__msa_selector.msa_file_selected.connect(self.__on_msa_selected)
        self.__ui.pruneButton.clicked.connect(self.__prune)
        self.__msa_viewer = MsaViewer()

        msa_viewer_layout = self.__ui.msaViewerContainerWidget.layout()

        assert msa_viewer_layout is not None, "The UI file has no layout for the msa viewer"

        msa_viewer_layout.replaceWidget(
            self.__ui.msaViewerWidget,
            self.__msa_viewer
        )
        self.__msa_viewer.show()

        self.__msa_viewer.range_selected.connect(self.__on_range_selected)

        self.__ui.saveResultsButtons.clicked.connect(self.__save_alignment)

    @pyqtSlot()
    def __on_range_selected(self):

        scope = self.__msa_viewer.selected_range

        if scope is None:
            return

        start, end = scope

        for _, cleaner in self.__cleaners:
            cleaner.set_scope(start, end)

    @pyqtSlot(QModelIndex, QModelIndex, 'QVector<int>')
    def __on_sequence_score_clicked(self, start : QModelIndex, stop: QModelIndex, roles : List[Qt.ItemDataRole] = []):
        if self.__scores_model is None or Qt.ItemDataRole.CheckStateRole not in roles:
            return

        self.__mask_msa_sequences()
        self.__update_stats()

    def __list_unwanted(self) -> Iterable[int]:

        assert(self.__scores_model)

        for score in self.__scores_model.scores:

            if score is None:
                continue
            for i,_ in enumerate(self.__scores_model.alignment):
                if not self.__scores_model.is_included(i):
                    yield i

    def __mask_msa_sequences(self):
        self.__msa_viewer.mask_sequences(self.__list_unwanted())

    @pyqtSlot(name="__on_score_changed")
    @throttle(1000)
    def __on_score_changed(self):
        if self.__scores_model is None:
            return

        self.__scores_model.update_scores(
            [
                ScoreEntry.from_result(name, result)
                for name, cleaner in self.__cleaners
                for result in [cleaner.score] if result is not None
            ]
        )

        self.__mask_msa_sequences()
        self.__update_stats()

    def __update_table(self, model : SequenceScoresModel) -> None:

        if self.__scores_model is not None:
            self.__scores_model.dataChanged.disconnect(self.__on_sequence_score_clicked)

        self.__scores_model = model
        self.__ui.scoresTable.setModel(model)
        self.__ui.scoresTable.resizeColumnsToContents()
        model.dataChanged.connect(self.__on_sequence_score_clicked)

    @pyqtSlot()
    def __prune(self):
        self.__update_alignment(self.__msa_viewer.get_new_alignment())
        show_info(self, "Alignment has been pruned")

    @pyqtSlot(name="__on_msa_selected")
    @with_error_handler()
    def __on_msa_selected(self):
        alignment = self.__msa_selector.alignment
        if alignment is None:
            raise Exception("The selected alignment could not be opened.")

        self.__update_alignment(alignment)

    def __update_alignment(self, alignment: MultipleSeqAlignment):

        context = ScoreContext.from_alignment(alignment)
        model = self.__scores_model
        scores = [
            ScoreEntry.from_result(name, result)
            for name, cleaner in self.__cleaners
            for result in [cleaner.score_alignment(context)] 
        ]

        if model is not None:
            prev_index = model.get_index_by_name()
            prev_meta = model.meta
            overrides = [
                index is not None and prev_meta[index]
                for seq in alignment
                for index in [prev_index.get(seq.id)]
            ]
            meta_override = lambda i,meta: overrides[i] or meta
        else:
            meta_override = None

        self.__update_table(
            SequenceScoresModel(
                alignment,
                scores,
                meta_override
            )
        )

        self.__msa_viewer.set_alignment(alignment)
        self.__mask_msa_sequences()
        self.__update_stats()

    @pyqtSlot(name = "__save_alignment")
    @with_error_handler()
    def __save_alignment(self):
        result_file,_ = QFileDialog.getSaveFileName(
            self,
            "Save File",
            "",
            "Multiple Sequence Alignment (*.fasta *.clustal)"
        )

        save_msa(result_file, self.__msa_viewer.get_new_alignment())

        show_info(
            self,
            "Saving Succeeded",
            f"Location: {result_file}"
        )

    def __update_stats(self) -> None:

        scores = self.__scores_model
        if scores is None:
            model_dict = {}
        else:
            total_seqs = len(scores.alignment)
            remaining_seqs = self.__msa_viewer.masked_alignment_seqs
            remaining_perc = round(100*remaining_seqs/total_seqs, 0)
            deleted_seqs = total_seqs - remaining_seqs
            deleted_seqs_perc = round(100*deleted_seqs/total_seqs, 0)
            alignment_length = scores.alignment.get_alignment_length()
            cleaned_length = self.__msa_viewer.masked_alignment_length
            cleaned_length_perc = round(100*cleaned_length / alignment_length, 0)
            deleted_pos = alignment_length - cleaned_length
            if deleted_seqs < 1:
                efficiency = 0
            else:
                efficiency = round(deleted_pos/deleted_seqs, 2)
            model_dict = {
                'Original Sequences': str(total_seqs),
                'Remaining Sequences': f"{remaining_seqs} ({remaining_perc}%)",
                'Deleted Sequences': f"{deleted_seqs} ({deleted_seqs_perc}%)",
                'Efficiency': f"{efficiency}",
                'Original Length': str(alignment_length),
                'Cleaned Length': f"{cleaned_length} ({cleaned_length_perc}%)"
            }

        self.__ui.statsTable.setModel(DictionaryModel(model_dict))
