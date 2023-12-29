from PyQt5.QtCore import QModelIndex, pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import QWidget
from typing import NamedTuple, Optional

from ...core.Qt.visual.NamedTupleEditor import MetaFieldOverrides, MetaFieldOverridesDict, namedtuple_eidtor
from ..data import Point2d, Series
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
    show_excluded : bool
    slope : float
    b : float

    def apply(self, x: float) -> float:
        return x*self.slope + self.b

    @staticmethod
    def update_values(
        series: Series,
        current: Optional['LinearFitValues'] = None
    ):

        points = series.values

        if current is None:
            x_min = min(point.x for point in points)
            x_max = max(point.x for point in points)
            show_excluded = True
        else:
            x_min = current.x_min
            x_max = current.x_max
            show_excluded = current.show_excluded

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
            show_excluded=show_excluded,
            slope = float(slope),
            b = float(b)
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
            ),
            'show_excluded': MetaFieldOverrides(
                display="Show excluded"
            )
        }

class LinearFitWidget(QWidget):

    lse_model_changed = pyqtSignal()

    def __init__(self, fit_meta: Optional[LinearFitMeta] = None):
        super().__init__()
        self.__series = None
        self.__model = None
        self.__lse_model = None
        self.__fit_meta = fit_meta or LinearFitMeta()
        self.__overrides = linear_fit_values_overrides(self.__fit_meta)

        self.__ui = Ui_LinearFitWidget()
        self.__ui.setupUi(self)

        self.__plot = Plot()

        self.layout().replaceWidget(
            self.__ui.plotWidget,
            self.__plot
        )

    def lse_model(self) -> Optional[LinearFitValues]:
        return self.__lse_model

    def set_x_min(self, x_min: float):
        model = self.__model

        assert model is not None, "Series is not set"
        
        value = model[0]

        assert value is not None, "Series is not set"

        model[0] = value._replace(x_min = x_min)

    @pyqtSlot(QModelIndex, QModelIndex)
    def __on_data_changed(self, start: QModelIndex, end: QModelIndex):

        assert self.__model is not None, "Bug in the code. Model should not be None"
        assert self.__series is not None, "Bug in the code. Series cannot be None"

        value = self.__model[0]

        if value == self.__lse_model:
            return

        self.__lse_model = self.__model[0] = LinearFitValues.update_values(
            self.__series,
            value
        )

        self.__plot_series()
        self.lse_model_changed.emit()

    def __plot_series(self):

        series = self.__series
        lse_model = self.__lse_model

        if series is None or lse_model is None:
            return self.__plot.set_series(SeriesSet())

        assert lse_model is not None, "There should always be one model"
        show_excluded = lse_model.show_excluded
        x_min = lse_model.x_min
        x_max = lse_model.x_max

        def included(point: Point2d) -> bool:
            return show_excluded \
                or (
                    point.x >= x_min and \
                    point.x <= x_max
                )

        line = [
            Point2d(
                x = point.x,
                y = lse_model.apply(point.x)
            )
            for point in series.values
                if included(point)
        ]

        self.__plot.set_series(
            SeriesSet(
                series=[
                    series.filter(included),
                    Series(
                        metadata=PlotMeta("LSE"),
                        values=line
                    )
                ],
                x_label=self.__fit_meta.x_axis_name,
                y_label=self.__fit_meta.y_axis_name
            )
        )

    def set_series(self, series : 'Series[PlotMeta]'):
        self.__series = series

        if self.__model is not None:
            self.__model.dataChanged.disconnect(
                self.__on_data_changed
            )

        self.__lse_model = LinearFitValues.update_values(series)
        self.__model = namedtuple_eidtor(
            self.__ui.statsTable,
            self.__lse_model,
            tuple_field_overrides=self.__overrides
        )

        self.__model.dataChanged.connect(self.__on_data_changed)

        self.__plot_series()
        self.lse_model_changed.emit()