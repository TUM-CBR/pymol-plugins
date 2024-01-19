from PyQt5.QtCore import pyqtSlot, QAbstractTableModel, QModelIndex, Qt
from PyQt5.QtWidgets import QWidget
from typing import Any, NamedTuple

from .Plot import Plot, SeriesSet
from .SeriesModel import SeriesModel
from .Ui_PlotWidget import Ui_PlotWidget

class SelectionState(NamedTuple):
    is_selected : bool = True

    def checked_state(self):
        if self.is_selected:
            return Qt.CheckState.Checked
        else:
            return Qt.CheckState.Unchecked
        
    @staticmethod
    def build_default(items: SeriesSet):
        num_series = len(items.series)
        return [
            SelectionState()
            for _i in range(0, num_series)
        ]

class SelectedSeriesModel(QAbstractTableModel):

    def __init__(
        self,
        series_set: SeriesSet
    ) -> None:
        super().__init__()

        self.__series_set = series_set
        self.__selected_state = SelectionState.build_default(series_set)

    def update_series_set(self, series_set: SeriesSet):
        old_length = len(self.__series_set.series)
        self.__series_set = series_set
        if len(series_set.series) != old_length:
            self.__selected_state = SelectionState.build_default(series_set)

        self.modelReset.emit()

    def rowCount(self, parent: Any = None) -> int:
        return len(self.__series_set.series)
    
    def columnCount(self, parent: Any = None) -> int:
        return 1
    
    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Vertical:
            return self.__series_set.series[section].metadata.name
        
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            assert section == 0, "This view model should have 1 column."
            return "Active"

        return super().headerData(section, orientation, role)
    
    def flags(self, index: QModelIndex) -> Qt.ItemFlags:
        return super().flags(index) | Qt.ItemFlag.ItemIsUserCheckable
    
    def __get_item(self, index: QModelIndex) -> SelectionState:
        return self.__selected_state[index.row()]
    
    def __set_item(self, index: QModelIndex, value: SelectionState):
        self.__selected_state[index.row()] = value
    
    def setData(self, index: QModelIndex, value: Any, role: int = Qt.ItemDataRole.EditRole) -> bool:

        if role == Qt.ItemDataRole.CheckStateRole:
            self.__set_item(
                index,
                self.__get_item(index)._replace(is_selected = value == Qt.CheckState.Checked)
            )
            return True

        return super().setData(index, value, role)
    
    def set_all(self, is_selected: bool):
        self.__selected_state = [
            value._replace(is_selected = is_selected)
            for value in self.__selected_state
        ]
        self.dataChanged.emit(
            self.index(0,0),
            self.index(len(self.__selected_state), 0)
        )
    
    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if not index.isValid():
            return None
        
        if role == Qt.ItemDataRole.CheckStateRole:
            return self.__get_item(index).checked_state()

        return None

    def current_selection(self) -> SeriesSet:
        return self.__series_set._replace(
            series = [
                series
                for selected, series in zip(self.__selected_state, self.__series_set.series) \
                    if selected.is_selected
            ]
        )

class PlotWidget(QWidget):

    def __init__(self) -> None:
        super().__init__()

        self.__ui = Ui_PlotWidget()
        self.__ui.setupUi(self)
        self.__x_scale = 1
        self.__y_scale = 1

        self.__plot = Plot(
            x_scale=self.__x_scale,
            y_scale=self.__y_scale
        )
        self.__selection_model = None

        plot_layout = self.__ui.plotContainer.layout()

        assert plot_layout is not None, "Bug in the code. No layout for plot container."

        plot_layout.replaceWidget(
            self.__ui.plotWidget,
            self.__plot
        )

        self.__ui.updateButton.clicked.connect(self.__on_update_clicked)
        self.__ui.allButton.clicked.connect(self.__set_all)
        self.__ui.noneButton.clicked.connect(self.__set_none)

    def set_scaling(self, x_scale: float = 1, y_scale: float = 1):
        self.__x_scale = x_scale
        self.__y_scale = y_scale
        self.__plot.set_scaling(x_scale, y_scale)

        model = self.__ui.dataTable.model()
        if model is None:
            return
        
        assert isinstance(model, SeriesModel), "The dataTable should have a SeriesModel"
        model.set_scaling(x_scale, y_scale)

    def __set_all(self):
        model = self.__selection_model

        if model is None:
            return

        model.set_all(True)

    def __set_none(self):
        model = self.__selection_model

        if model is None:
            return

        model.set_all(False)

    def set_series(self, series: SeriesSet):

        model = self.__selection_model
        if model is None:
            self.__selection_model = SelectedSeriesModel(series)
        else:
            model.update_series_set(series)

        self.__ui.selectedSeriesTable.setModel(self.__selection_model)
        self.__update()

    @pyqtSlot()
    def __on_update_clicked(self):
        self.__update()

    def __update(self):

        selection_model = self.__selection_model

        if selection_model is None:
            return

        series = selection_model.current_selection()
        self.__plot.set_series(series)
        self.__ui.dataTable.setModel(
            SeriesModel(
                series.series,
                render_series_header=lambda i,s: s.metadata.name,
                x_scale=self.__x_scale,
                y_scale=self.__y_scale
            )
        )

    