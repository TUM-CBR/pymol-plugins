from PyQt5.QtWidgets import QWidget
from typing import Iterable, Tuple

from ..data import *
from ..kinetics import *
from .LinearFitWidget import LinearFitMeta, LinearFitWidget
from .Plot import Plot, PlotMeta, Point2d, Series, SeriesSet
from .Ui_KineticsOptimizer import Ui_KineticsOptimizer


# km should be the concentration at which we have half v_max
# v_max := 1/v vs 1/s, km = v_max * m, v_max = 1/y-intercept
def velocity(v_max: float, beta: float, k_si: float, k_m: float, conc: float):

    num = (1 + conc/k_si * beta) * conc
    den = (1 + conc/k_si) * conc + k_m

    return v_max * num / den

def approximate_vmax_and_km(runs: List[List[Point2d]]) -> Tuple[float, float]:
    xs = [1/max(point.x for point in series) for series in runs]
    ys = [1/avg([point.y for point in series]) for series in runs]
    m,b = lse(xs, ys)
    return (1/b, m/b)

def aproximate_k_si_and_beta(
    runs : List[List[Point2d]],
    v_max: float
    ) -> Tuple[float, float]:

    scaled_runs = [
        Point2d(
            x = max(point.x for point in series),
            y = avg([point.y for point in series])
        )
        for series in runs
    ]

    min_conc = min(
        (point.x for point in scaled_runs if point.y < v_max),
    )

    valid_runs = [point for point in scaled_runs if point.x >= min_conc and point.y < v_max]

    x_values = [
        1 / point.x
        for point in valid_runs
    ]

    y_values = [
        v/(v_max - v)
        for point in valid_runs
        for v in [point.y]
    ]

    m, b = lse(x_values, y_values)
    return (b, m)

def vel_vs_conc_theory(runs: List[List[Point2d]]) -> List[Point2d]:
    runs.sort(key=lambda run: max(point.x for point in run))
    low_conc = max(2, int(len(runs) / 4))
    runs_low = runs[0:low_conc]
    v_max, k_m = approximate_vmax_and_km(runs_low)

    k_si, beta = aproximate_k_si_and_beta(runs, v_max)

    return [
        Point2d(
            x = point.x,
            y = velocity(v_max, beta, k_si, k_m, point.x)
        )
        for series in runs
        for point in series
    ]


def as_vmax_vs_km(runs : 'List[Series[RunMetadata]]') -> 'Series[PlotMeta]':
    return Series(
        metadata=PlotMeta(name = "data"),
        values = [
            Point2d(
                x = 1/series.metadata.concentration,
                y = 1/avg(point.y for point in series.values)
            )

            for series in runs
        ]
    )

K_TAB_1 = "Vmax and Km"
K_TAB_2 = "Beta and Ki"

class KineticsOptimizer(QWidget):

    def __init__(
        self,
        runs: KineticsRuns
    ) -> None:
        super().__init__()

        self.__ui = Ui_KineticsOptimizer()
        self.__ui.setupUi(self)

        self.__plot_widget = Plot()
        self.layout().replaceWidget(
            self.__ui.plotWidget,
            self.__plot_widget
        )

        self.__vmax_widget = LinearFitWidget(
            LinearFitMeta(
                x_axis_name="1/[S]",
                y_axis_name="1/v"
            )
        )
        self.__ui.initalValuesWidget.insertTab(
            0,
            self.__vmax_widget,
            K_TAB_1
        )

        self.__beta_widget = LinearFitWidget()

        self.__ui.initalValuesWidget.insertTab(
            1,
            self.__beta_widget,
            K_TAB_2
        )

        self.__runs = runs
        self.__velocity_vs_conc = as_vel_vs_conc_series(runs)

        self.__render_plots()

    def __render_vmax_vs_km(self):

        values = [
            Point2d(
                x = 1/point.x,
                y = 1/point.y
            )
            for series in self.__velocity_vs_conc
            for point in series.values
        ]

        values.sort(key=lambda point: point.x)
        self.__vmax_widget.set_series(
            Series(
                PlotMeta(name=K_TAB_1),
                values
            )
        )

    def __render_plots(self):

        self.__render_vmax_vs_km()
        #exp = as_vel_vs_conc_series(self.__runs)
        #self.__plot_widget.set_series(exp)