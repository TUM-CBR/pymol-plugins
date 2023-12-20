from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QVBoxLayout, QWidget 
from typing import List, NamedTuple

from ..data import *

class PlotMeta(NamedTuple):
    name : str

class SeriesSet(NamedTuple):
    series : List['Series[PlotMeta]'] = []
    x_label : str = "x"
    y_label : str = "y"

class Plot(QWidget):

    def __init__(self):
        super().__init__()
        self.__series = SeriesSet()

    def set_series(self, series: SeriesSet):
        self.__series = series
        self.__render_series()

    def __render_series(self):
        figure = Figure()
        canvas = FigureCanvas(figure)
        ax = figure.add_subplot(111)
        series_set = self.__series

        for series in series_set.series:
            ax.plot(
                [point.x for point in series.values],
                [point.y for point in series.values],
                label=series.metadata.name
            )

        ax.set_xlabel(series_set.x_label)
        ax.set_ylabel(series_set.y_label)
        ax.legend()
        layout = QVBoxLayout(self)
        layout.addWidget(canvas)

        self.setLayout(layout)