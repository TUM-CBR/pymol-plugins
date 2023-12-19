from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QVBoxLayout, QWidget 
from typing import List, NamedTuple

class Point2d(NamedTuple):
    x : float
    y : float

class Series(NamedTuple):
    name : str
    values : List[Point2d]

class Plot(QWidget):

    def __init__(self):
        super().__init__()
        self.__series = []

    def set_series(self, series: List[Series]):
        self.__series = series
        self.__render_series()

    def __render_series(self):
        figure = Figure()
        canvas = FigureCanvas(figure)
        ax = figure.add_subplot(111)

        for series in self.__series:
            ax.plot(
                [point.x for point in series.values],
                [point.y for point in series.values],
                label=series.name
            )

        ax.set_xlabel('time')
        ax.set_ylabel('abs')
        ax.legend()
        layout = QVBoxLayout(self)
        layout.addWidget(canvas)

        self.setLayout(layout)