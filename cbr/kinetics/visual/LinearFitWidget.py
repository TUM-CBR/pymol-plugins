from PyQt5.QtCore import QModelIndex, pyqtSlot
from PyQt5.QtWidgets import QWidget
from typing import NamedTuple, Optional

from ...core.Qt.visual.NamedTupleEditor import MetaFieldOverrides, MetaFieldOverridesDict, namedtuple_eidtor
from ..data import Series
from ..math import *

from .Plot import Plot, PlotMeta, SeriesSet
from .Ui_LinearFitWidget import Ui_LinearFitWidget

class LinearFitMeta(NamedTuple):
    slope_name : str = "m"
    intercept_name : str = "b"
    x_axis_name : str = "x"
    y_axis_name : str = "y"

def lse_series(
    series: Series,
    x_min: float,
    x_max: float
):
    points = [
        point
        for point in series.values if point.x >= x_min and point.x <= x_max
    ]
    return lse(
        [point.x for point in points],
        [point.y for point in points]
    )

class LinearFitValues(NamedTuple):

    x_min : float
    x_max : float
    slope : float
    b : float

    @staticmethod
    def update_values(
        series: Series,
        current: Optional['LinearFitValues'] = None
    ):

        points = series.values

        if current is None:
            x_min = min(point.x for point in points)
            x_max = max(point.x for point in points)
        else:
            x_min = current.x_min
            x_max = current.x_max

        points = [
            point
            for point in points
                if point.x >= x_min and point.x <= x_max
        ]
        slope, b = lse(
            [point.x for point in points],
            [point.y for point in points]
        )

        return LinearFitValues(
            x_min = x_min,
            x_max = x_max,
            slope = slope,
            b = b
        )

def linear_fit_values_overrides(meta: LinearFitMeta) -> MetaFieldOverridesDict:

    return {
            'x_min': MetaFieldOverrides(
                display=f"min ({meta.x_axis_name})"
            ),
            'x_max': MetaFieldOverrides(
                display=f"max ({meta.x_axis_name})"
            ),
            'slope': MetaFieldOverrides(
                display=meta.slope_name,
                readonly=True
            ),
            'b': MetaFieldOverrides(
                display=meta.intercept_name,
                readonly=True
            )
        }

class LinearFitWidget(QWidget):

    def __init__(self, fit_meta: Optional[LinearFitMeta] = None):
        super().__init__()
        self.__series = None
        self.__model = None
        self.__value = None
        self.__fit_meta = fit_meta or LinearFitMeta()
        self.__overrides = linear_fit_values_overrides(self.__fit_meta)

        self.__ui = Ui_LinearFitWidget()
        self.__ui.setupUi(self)

        self.__plot = Plot()

        self.layout().replaceWidget(
            self.__ui.plotWidget,
            self.__plot
        )

    @pyqtSlot(QModelIndex, QModelIndex)
    def __on_data_changed(self, start: QModelIndex, end: QModelIndex):

        assert self.__model is not None, "Bug in the code. Model should not be None"
        assert self.__series is not None, "Bug in the code. Series cannot be None"

        value = self.__model[0]

        if value == self.__value:
            return

        self.__value = self.__model[0] = LinearFitValues.update_values(
            self.__series,
            value
        )

    def set_series(self, series : 'Series[PlotMeta]'):
        self.__series = series

        if self.__model is not None:
            self.__model.dataChanged.disconnect(
                self.__on_data_changed
            )

        self.__model = namedtuple_eidtor(
            self.__ui.statsTable,
            LinearFitValues.update_values(
                series
            ),
            tuple_field_overrides=self.__overrides
        )

        self.__plot.set_series(
            SeriesSet(
                series=[series],
                x_label=self.__fit_meta.x_axis_name,
                y_label=self.__fit_meta.y_axis_name
            )
        )

        self.__model.dataChanged.connect(self.__on_data_changed)