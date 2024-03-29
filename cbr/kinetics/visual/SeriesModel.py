from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt
from typing import Any, Callable, Generic, List

from ..data import Series, TMeta

def default_series_header(i : int, series: 'Series[Any]') -> str:
    return f"Series {i}"

class SeriesModel(QAbstractTableModel, Generic[TMeta]):

    def __init__(
        self,
        series : List['Series[TMeta]'],
        x_scale : float,
        y_scale : float,
        render_series_header : Callable[[int, 'Series[TMeta]'], str] = default_series_header
    ) -> None:
        super().__init__()
        self.__render_series_header = render_series_header
        self.__x_scale = x_scale
        self.__y_scale = y_scale

        self.__set_series(series)

    def __set_series(self, series: List['Series[TMeta]']):
        self.__series = series
        self.__row_count = max(len(series.values) for series in series)
        self.__column_headers = [
            self.__render_series_header(i, series)
            for i,series in enumerate(series)
        ]

    def rowCount(self, parent: Any = None) -> int:
        return self.__row_count

    def columnCount(self, parent: Any = None) -> int:
        return len(self.__column_headers)
    
    def set_scaling(self, x_scale: float, y_scale: float):
        self.__x_scale = x_scale
        self.__y_scale = y_scale
        self.modelReset.emit()

    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:

        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            return self.__column_headers[section]

        return super().headerData(section, orientation, role)

    def __display_role_data(self, index: QModelIndex) -> Any:

        return f"{round(self.__series[index.column()].values[index.row()].y * self.__y_scale, 8)}"

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if not index.isValid():
            return None

        if role == Qt.ItemDataRole.DisplayRole:
            return self.__display_role_data(index)