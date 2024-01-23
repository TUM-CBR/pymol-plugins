from PyQt5.QtCore import pyqtSlot, QAbstractTableModel, QModelIndex, QObject, Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QSpinBox, QTableView
from typing import Any, Dict, Iterable, List, Optional

from cbr.core.Qt.visual.JsonRecordsTable import Record

from ...core import color
from ...core import dictutils
from ...core import optional
from ...core.Qt.visual.JsonRecordsTable import JsonRecordsModel

Fragments = Dict[str, List[str]]

colors = [
    QColor(*rgb)
    for rgb in color.generate_colors(10)
]

def get_shuffling_points(seq: str, count: int) -> List[int]:
    seg_size = len(seq) / (count + 1)
    return [
        int(seg_size*i)
        for i in range(1, count + 1)
    ]

def get_all_shuffling_points(
    seqs: Dict[str, str],
    count : int
) -> Iterable[Dict[str, int]]:

    sps = dict(
        (name, get_shuffling_points(seq, count))
        for name,seq in seqs.items()
    )

    for i in range(0, count):
        yield dict(
            (name, sps[name][i])
            for name in seqs.keys()
        )

def get_segments(seq: str, shuffling_points: List[int]) -> Iterable[str]:
    start = 0

    for ix in shuffling_points:
        yield seq[start:ix]
        start = ix

def get_lengths(seq_segments: Dict[str, List[str]]) -> List[int]:

    sizes = [
        [len(seg) for seg in segments]
        for segments in seq_segments.values()
    ]
    count = len(sizes[0])
    assert all(len(size) == count for size in sizes), "All sequences must have the same number of shuffling points"

    return [
        max(size[i] for size in sizes)
        for i in range(0, count)
    ]

def add_padding(seq: str, target: int) -> str:
    ammount = max(target - len(seq), 0)
    return " "*ammount + seq

def as_positions(lengths: List[int]) -> Iterable[int]:
    offset = 0
    for length in lengths:
        yield offset + length
        offset += length

class SequencesBuilderModel(QAbstractTableModel):

    def __init__(
        self,
        sequences: Dict[str, str]
    ):
        super().__init__()
        self.__seqs_names = list(sequences.keys())
        self.__sequences = sequences
        self.__columns = 0
        self.__fragments = None

    def sequences(self) -> Dict[str, str]:
        return self.__sequences
    
    def fragments(self) -> Optional[Fragments]:
        return self.__fragments

    def set_shuffling_points(self, shuffling_points: Dict[str, List[int]]):

        self.__fragments = segments = dict(
            (name, list(get_segments(self.__sequences[name], points)))
            for name, points in shuffling_points.items()
        )

        desired_lengths = get_lengths(segments)
        self.__columns = sum(desired_lengths)

        self.__display_sequences = dict(
            (
                name,
                "".join(add_padding(seq, desired_lengths[i]) for i,seq in enumerate(seqs))
            )
            for name, seqs in segments.items()
        )

        self.__positions = list(as_positions(desired_lengths))
        self.modelReset.emit()

    def rowCount(self, parent: Any = None) -> int:
        return len(self.__sequences)
    
    def columnCount(self, parent: Any = None) -> int:
        return self.__columns
    
    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if orientation == Qt.Orientation.Vertical and role == Qt.ItemDataRole.DisplayRole:
            return self.__seqs_names[section]
        return super().headerData(section, orientation, role)
    
    def __get_value(self, index: QModelIndex) -> str:
        seq_name = self.__seqs_names[index.row()]
        seq = self.__display_sequences[seq_name]

        return seq[index.column()]
    
    def __get_color(self, index: QModelIndex) -> QColor:

        color_index = next(
            i
            for i,position in enumerate(self.__positions) if index.column() < position
        )

        return colors[color_index % len(colors)]
    
    def data(
        self,
        index: QModelIndex,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if not index.isValid():
            return
        
        if role == Qt.ItemDataRole.DisplayRole:
            return self.__get_value(index)
        elif role == Qt.ItemDataRole.BackgroundColorRole:
            return self.__get_color(index)

        return None


class ShufflingPointsModel(JsonRecordsModel):

    def __init__(self, data: List[Record]):
        super().__init__(data, False)

    def get_shuffling_points(self) -> Dict[str, List[int]]:

        return dictutils.merge(
            lambda vs: list(optional.assert_values(vs)),
            *self.records()
        )


class SequencesPositionEditor(QObject):

    def __init__(
        self,
        points_table: QTableView,
        sequences_table: QTableView,
        shuffling_points_box: QSpinBox
    ) -> None:
        super().__init__()
        self.__points_table = points_table
        self.__sequences_table = sequences_table
        self.__shuffling_points_box = shuffling_points_box
        self.__shuffling_points_model = ShufflingPointsModel([])
        self.__shuffling_points_model.dataChanged.connect(self.__on_shuffling_points_changed)
        self.__points_table.setModel(self.__shuffling_points_model)
        self.__sequences_builder_model = None
        self.__shuffling_points_box.valueChanged.connect(self.__on_shuffling_points_count_changed)

    @pyqtSlot(QModelIndex, QModelIndex)
    def __on_shuffling_points_changed(self, start: Any, end: Any):
        self.__update_shuffling_points()

    def __update_shuffling_points(self):
        seqs_model = self.__sequences_builder_model
        shuffling_points_model = self.__shuffling_points_model

        if seqs_model is None:
            return
        
        seqs_model.set_shuffling_points(shuffling_points_model.get_shuffling_points())
        self.__sequences_table.resizeColumnsToContents()

    @pyqtSlot()
    def __on_shuffling_points_count_changed(self):

        seqs_model = self.__sequences_builder_model

        if seqs_model is None:
            return

        self.__set_shuffling_points_count(seqs_model.sequences())
        self.__update_shuffling_points()

    def __set_shuffling_points_count(self, seqs: Dict[str, str]):
        count = self.__shuffling_points_box.value()
        shuffling_points = list(get_all_shuffling_points(seqs, count))
        self.__shuffling_points_model.set_data(shuffling_points)


    def set_sequences(self, seqs: Dict[str, str]):

        self.__sequences_builder_model = model = SequencesBuilderModel(seqs)
        self.__sequences_table.setModel(model)
        self.__set_shuffling_points_count(seqs)

        model.set_shuffling_points(self.__shuffling_points_model.get_shuffling_points())
        self.__sequences_table.resizeColumnsToContents()

    def fragments(self) -> Optional[Fragments]:
        seqs_model = self.__sequences_builder_model

        if seqs_model is None:
            return None
        
        return seqs_model.fragments()