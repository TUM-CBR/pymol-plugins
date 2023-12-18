from PyQt5.QtGui import QColor, QPaintEvent, QPen, QPixmap
from PyQt5.QtWidgets import QWidget 
from typing import Callable, List, NamedTuple, Optional, Tuple

from ...core import color
from ...core.Qt.QtGui import paint
from ...control import pairwise

class Point2d(NamedTuple):
    x : float
    y : float

class Series(NamedTuple):
    name : str
    values : List[Point2d]

class Plot(QWidget):

    def __init__(self):
        self.__series = []
        self.__x_domain = None
        self.__y_domain = None
        self.__plot_res = 1000
        self.__pen_width = 5
        self.__pixmap = None

    def set_series(self, series: List[Series]):
        self.__series = series
        self.__render_series()

    def __domain_of(
        self,
        key: Callable[[Point2d], float],
        override : Optional[Tuple[float, float]]
    ) -> Tuple[float, float]:
        if override is not None:
            return override
        elif len(self.__series) == 0:
            return (0,0)

        min_x = min(
            key(point)
            for series in self.__series
            for point in series.values
        )

        max_x = max(
            key(point)
            for series in self.__series
            for point in series.values
        )

        return (min_x, max_x)

    @property
    def x_domain(self):
        return self.__domain_of(lambda p: p.x, self.__x_domain)

    @x_domain.setter
    def x_domain(self, value: Tuple[float, float]):
        self.__x_domain = value

    @property
    def y_domain(self):
        return self.__domain_of(lambda p: p.y, self.__y_domain)

    @y_domain.setter
    def y_domain(self, value: Tuple[float, float]):
        self.__y_domain = value

    def get_color_for(self, i : int) -> QColor:
        return QColor(color.get_qt_color(i))

    def __render_series(self):

        self.__pixmap = pixmap = QPixmap(self.__plot_res, self.__plot_res)

        with paint(pixmap) as manager:

            painter = manager.painter
            x_min, x_max = self.x_domain
            x_size = x_max - x_min
            x_res = self.__plot_res
            x_scale = abs(x_res / x_size)
            def x_mapper(x: float):
                x = int((x - x_min) * x_scale)

                if x >= x_res:
                    return x_res - 1
                elif x < 0:
                    return 0
                else:
                    return x
            
            y_min, y_max = self.y_domain
            y_size = y_max - y_min
            y_res = self.__plot_res
            y_scale = abs(x_res / y_size)
            def y_mapper(y: float):
                y = int((y - y_min) * y_scale)

                if y >= y_res:
                    return y_res - 1
                elif y < 0:
                    return 0
                else:
                    return y

            for i,series in enumerate(self.__series):
                pen = QPen(self.get_color_for(i))
                pen.setWidth(self.__pen_width)

                #itertools.pairwise not in python 3.7
                for prev, current in pairwise(series.values):
                    painter.drawLine(
                        x_mapper(prev.x),
                        y_mapper(prev.y),
                        x_mapper(current.x),
                        y_mapper(current.y)
                    )

    def paintEvent(self, a0: QPaintEvent) -> None:

        super().paintEvent(a0)

        pixmap = self.__pixmap

        if pixmap is None:
            return

        with paint(self) as manager:
            manager.painter.drawPixmap(
                0,
                0,
                self.width(),
                self.height(),
                pixmap,
                0,
                0,
                pixmap.width(),
                pixmap.height()
            )