from pymol import cmd
from PyQt5.QtCore import QAbstractTableModel, QItemSelection, QModelIndex, QObject, QThread, Qt, pyqtSignal, pyqtSlot
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QWidget
from typing import Any, cast, List, Optional, Tuple

from ...core.Context import Context
from ...core.pymol.algorithms import SpatialAlignment, align_by_rmsd
from ...core.pymol.structure import StructureSelection
from ...core.Qt.QtWidgets import with_error_handler
from ...core.visual import as_structure_selector
from ...support.display.sequence import RESIDUE_COLORS
from .Ui_FrankenProt import Ui_FrankenProt

class AlignRmsdWorker(QObject):

    rmsd_started = pyqtSignal()

    rmsd_completed = pyqtSignal(object, object)

    @pyqtSlot(object, object)
    def run_rmsd_align(self, s1: StructureSelection, s2: StructureSelection):

        self.rmsd_started.emit()

        try:
            result = align_by_rmsd(s1, s2)
            self.rmsd_completed.emit(result, None)
        except Exception as e:
            self.rmsd_completed.emit(None, e)


class FrankenProtModel(QAbstractTableModel):

    BASE_SEQ_ROW_IX = 0
    FRAGMENT_SEQ_ROW_IX = 1
    DISTANCE_ROW_IX = 2

    ROWS = [
        BASE_SEQ_ROW_IX,
        FRAGMENT_SEQ_ROW_IX,
        DISTANCE_ROW_IX
    ]

    DISTANCE_CAP = 5

    def __init__(self, alignment: SpatialAlignment):
        super().__init__()
        self.__alignment = alignment
        self.__distances = [
            self.__get_distance_for_position(i)
            for i in range(alignment.alignment.get_alignment_length())
        ]

    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.ROWS)
    
    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        return self.__alignment.alignment.get_alignment_length()
    
    def __row_header(self, ix: int) -> Optional[str]:

        if ix == self.BASE_SEQ_ROW_IX:
            return cast(Any, self.__alignment.alignment[0]).id
        elif ix == self.FRAGMENT_SEQ_ROW_IX:
            return cast(Any, self.__alignment.alignment[1]).id
        elif ix == self.DISTANCE_ROW_IX:
            return "Distance Error"
        else:
            return None
        
    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if role != Qt.ItemDataRole.DisplayRole:
            return super().headerData(section, orientation, role)
        
        if orientation == Qt.Orientation.Vertical:
            return self.__row_header(section)
        
        return super().headerData(section, orientation, role)
    
    def __position_from_index(self, index: QModelIndex) -> int:
        return index.column()

    def __sequence_and_position(self, index: QModelIndex) -> Optional[Tuple[int, int]]:

        row = index.row()

        if row == self.BASE_SEQ_ROW_IX:
            return (0, self.__position_from_index(index))
        elif row == self.FRAGMENT_SEQ_ROW_IX:
            return (1, self.__position_from_index(index))
        else:
            return None
        
    def __get_residue(self, index: QModelIndex) -> Optional[str]:
        res_position = self.__sequence_and_position(index)

        if res_position is None:
            return None
        
        seq, pos = res_position
        return self.__alignment.alignment[seq]._seq[pos] # type: ignore
    
    def __get_distance_for_position(self, pos: int) -> Optional[float]:
        resv1 = self.__alignment.resv_maps[0].get(pos)
        resv2 = self.__alignment.resv_maps[1].get(pos)

        if resv1 is None or resv2 is None:
            return None
        
        return cast(Any, self.__alignment.distance_by_resv(resv1, resv2))

    def __get_distance(self, index: QModelIndex) -> Optional[float]:
        pos = self.__position_from_index(index)
        return self.__distances[pos]
    
    def __data_role(self, index: QModelIndex) -> Optional[str]:
        resv = self.__get_residue(index)
        if resv is not None:
            return resv
        
        distance = self.__get_distance(index)

        if distance is not None:
            return str(round(distance, 2))
        
    def __distance_color(self, value: float):

        value = min(value, self.DISTANCE_CAP)
        green = int(255 * (1 - value/self.DISTANCE_CAP))
        red = int(255 * value/self.DISTANCE_CAP)

        return QColor(red,green,0)
    
    def __color_role(self, index: QModelIndex) -> Optional[QColor]:

        residue = self.__get_residue(index)

        if residue is not None:
            return RESIDUE_COLORS.get(residue)
        
        distance = self.__get_distance(index)

        return self.__distance_color(distance if distance is not None else self.DISTANCE_CAP)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if not index.isValid():
            return None
        
        if role == Qt.ItemDataRole.DisplayRole:
            return self.__data_role(index)
        elif role == Qt.ItemDataRole.BackgroundColorRole:
            return self.__color_role(index)
        else:
            return None
        
    def on_columns_selected(self, indexes: List[QModelIndex]):

        if len(indexes) == 0:
            return

        columns = set(self.__position_from_index(i) for i in indexes)
        base_resv = [
            resv
            for c in columns
            for resv in [self.__alignment.resv_maps[0].get(c)] if resv is not None
        ]
        fragment_resv = [
            resv
            for c in columns
            for resv in [self.__alignment.resv_maps[1].get(c)] if resv is not None
        ]

        selections = self.__alignment.distance_matrix.selections
        base_model = selections[0]
        base_sel = base_model.selection

        fragment_model = selections[1]
        fragment_sel = fragment_model.selection

        cmd.show(representation='cartoon', selection=base_sel)
        cmd.color("blue", base_sel)

        if len(base_resv) > 0:
            cmd.hide(selection=base_model.residue_selection(base_resv))

        cmd.hide(selection=fragment_sel)

        if len(fragment_resv) > 0:
            start = min(fragment_resv)
            end = max(fragment_resv)
            cmd.show(representation='cartoon', selection=fragment_model.residue_selection([start - 1, end + 1] + fragment_resv))
            cmd.color("green", fragment_sel)

class FrankenProt(QWidget):

    start_rmsd = pyqtSignal(object, object)

    def __init__(self, context: Context):
        super().__init__()
        self.__ui = Ui_FrankenProt()
        self.__ui.setupUi(self)

        self.__base_structure = as_structure_selector(
            self.__ui.baseStructureCombo,
            self.__ui.refreshButton
        )

        self.__fragment_structure = as_structure_selector(
            self.__ui.fragmentsStructureCombo,
            self.__ui.refreshButton
        )

        self.__ui.selectButton.clicked.connect(self.__on_select_clicked)
        self.__ui.statusProgress.setVisible(False)

        self.__align_thread = QThread()
        self.__align_worker = AlignRmsdWorker()
        self.__align_worker.rmsd_started.connect(self.__on_rmsd_started)
        self.__align_worker.rmsd_completed.connect(self.__on_rmsd_completed)
        self.start_rmsd.connect(self.__align_worker.run_rmsd_align)
        self.__align_worker.moveToThread(self.__align_thread)
        self.__align_thread.start()

        self.__model: Optional[FrankenProtModel] = None

    @pyqtSlot(QItemSelection, QItemSelection)
    def __on_selection_changed(self, selected: QItemSelection, deselected: QItemSelection):

        model = self.__model

        if model is None:
            return

        model.on_columns_selected(self.__selection_model.selectedIndexes())

    @pyqtSlot()
    def __on_rmsd_started(self):
        self.__ui.selectButton.setEnabled(False)
        self.__ui.statusProgress.setVisible(True)

    @pyqtSlot(object, object, name="__on_rmsd_completed")
    @with_error_handler()
    def __on_rmsd_completed(
        self,
        result: Optional[SpatialAlignment],
        exn: Optional[Exception]
    ):
        self.__ui.selectButton.setEnabled(True)
        self.__ui.statusProgress.setVisible(False)

        if exn is not None:
            raise exn
        
        if result is not None:
            self.__model = FrankenProtModel(result)
            self.__ui.sequencesTable.setModel(self.__model)
            selection_model = self.__ui.sequencesTable.selectionModel()
            assert selection_model is not None, "Selection model unexpectedly None"
            selection_model.selectionChanged.connect(self.__on_selection_changed)
            self.__selection_model = selection_model
            self.__ui.sequencesTable.resizeColumnsToContents()

    @pyqtSlot(name="__on_select_clicked")
    @with_error_handler()
    def __on_select_clicked(self):

        base_structure = self.__base_structure.currentSelection
        fragment_structure = self.__fragment_structure.currentSelection

        if base_structure is None or fragment_structure is None:
            raise Exception("You must select both structures")
        
        cmd.align(
            fragment_structure.selection,
            base_structure.selection
        )

        self.start_rmsd.emit(base_structure, fragment_structure)