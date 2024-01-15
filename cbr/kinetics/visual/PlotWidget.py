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

class SelectedSeriesModel(QAbstractTableModel):

    def __init__(self, series_set: SeriesSet) -> None:
        super().__init__()

        self.__series_set = series_set
        self.__selected_state = [
            SelectionState()
            for _ in range(0, len(series_set.series))
        ]

    def rowCount(self, parent = None) -> int:
        return len(self.__series_set.series)
    
    def columnCount(self, parent = None) -> int:
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
    
    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if not index.isValid():
            return None
        
        if role == Qt.ItemDataRole.CheckStateRole:
            return self.__get_item(index).checked_state()

        return super().data(index, role)

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

        self.__plot = Plot()
        self.__selection_model = None

        plot_layout = self.__ui.plotContainer.layout()

        assert plot_layout is not None, "Bug in the code. No layout for plot container."

        plot_layout.replaceWidget(
            self.__ui.plotWidget,
            self.__plot
        )

        self.__ui.updateButton.clicked.connect(self.__on_update_clicked)

    def set_series(self, series: SeriesSet):
        self.__selection_model = SelectedSeriesModel(series)
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
                render_series_header=lambda i,s: s.metadata.name
            )
        )

    