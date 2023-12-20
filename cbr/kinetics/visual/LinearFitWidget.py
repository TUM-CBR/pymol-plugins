from PyQt5.QtWidgets import QWidget
from typing import Callable, NamedTuple, Optional

from ...core.Qt.visual.NamedTupleEditor import NamedTupleEditorModel
from ..data import Series
from ..math import *

from .Plot import Plot
from .Ui_LinearFitWidget import Ui_LinearFitWidget

class LinearFitMeta(NamedTuple):
    slope_name : str = "m"
    slop_mapping : Callable[[float], float] = lambda x: x
    intercept_name : str = "b"
    intercept_mapping : Callable[[float], float] = lambda b: b
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

class LinearFitModel(NamedTupleEditorModel):
    pass

class LinearFitWidget(QWidget):

    def __init__(self, fit_meta: Optional[LinearFitMeta] = None):
        super().__init__()
        self.__series = None

        self.__ui = Ui_LinearFitWidget()
        self.__ui.setupUi(self)

        self.__plot = Plot()

    def set_series(self, series : Series):
        self.__series = series