from PyQt5.QtGui import QHideEvent, QShowEvent
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QVBoxLayout, QWidget 
from typing import List, NamedTuple, Optional

from ..data import *

class PlotMeta(NamedTuple):
    name : str

class SeriesSet(NamedTuple):
    series : List['Series[PlotMeta]'] = []
    x_label : str = "x"
    y_label : str = "y"

class Plot(QWidget):

    def __init__(
        self,
        x_scale : float = 1,
        y_scale : float = 1
    ):
        super().__init__()
        self.__series = SeriesSet()
        self.__plot_canvas = FigureCanvas()
        self.__plot_layout = QVBoxLayout(self)
        layout = self.__plot_layout
        layout.addWidget(self.__plot_canvas)
        self.setLayout(layout)
        self.__x_scale = x_scale
        self.__y_scale = y_scale
        self.__should_render = False

    def set_series(self, series: SeriesSet):
        self.__series = series
        self.__render_series()

    def set_scaling(self, x_scale: float, y_scale: float):
        self.__x_scale = x_scale
        self.__y_scale = y_scale
        self.__render_series()

    def hideEvent(self, a0: Optional[QHideEvent]) -> None:
        self.__should_render = False
        return super().hideEvent(a0)
    
    def showEvent(self, a0: Optional[QShowEvent]) -> None:
        self.__should_render = True
        self.__render_series()
        return super().showEvent(a0)
    
    def __update_figure(self, figure: Figure):

        old_canvas = self.__plot_canvas
        self.__plot_canvas = FigureCanvas(figure=figure)
        self.__plot_layout.replaceWidget(
            old_canvas,
            self.__plot_canvas
        )

    def __render_series(self):

        if not self.__should_render:
            return

        figure = Figure()
        ax = figure.add_subplot(111)
        series_set = self.__series

        for series in series_set.series:
            ax.plot(
                [point.x * self.__x_scale for point in series.values],
                [point.y * self.__y_scale for point in series.values],
                label=series.metadata.name
            )

        ax.set_xlabel(series_set.x_label)
        ax.set_ylabel(series_set.y_label)
        ax.legend()

        self.__update_figure(figure)
