from PyQt5.QtWidgets import QWidget
from typing import Iterable, Tuple

from ..data import *
from .LinearFitWidget import LinearFitWidget
from .Plot import Plot, PlotMeta, Point2d, Series, SeriesSet
from .Ui_KineticsOptimizer import Ui_KineticsOptimizer

def is_baseline_run(run : KineticsRun) -> bool:
    return int(run.run_metadata.concentration * 1000) == 0

def run_to_vel_vs_conc(
    global_attributes : GlobalAttributes,
    run : KineticsRun,
    baseline : float = 0
) -> Iterable[Point2d]:
    conc_units = global_attributes.concentration_units
    intervals = global_attributes.measurement_interval
    conc = run.run_metadata.concentration
    molar_extinction = global_attributes.molar_extinction
    distance = global_attributes.distance
    prev_product_conc = 0

    for value in run.data:

        # Remove the baseline noise
        value -= baseline

        # abs = conc * molar_extinction * distance
        # => conc = abs / (molar_extinction * distance)
        product_conc = (value / (molar_extinction * distance)) / conc_units

        current_vel = abs(prev_product_conc - product_conc) / intervals
        prev_product_conc = product_conc

        # We assume that for every molecule of the product
        # we consume a molecule of the substrate
        yield Point2d(
            x = conc - product_conc,
            y = current_vel# / conc
        )

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

def as_vel_vs_conc_series(runs: KineticsRuns) -> List['Series[RunMetadata]']:

    baseline_run = next(
        (run.data for run in runs.runs if is_baseline_run(run)),
        [0.0]
    )
    baseline = avg(baseline_run)

    return [
        Series(
            metadata = run.run_metadata,
            values = list(run_to_vel_vs_conc(runs.global_attributes, run, baseline))
        )
        for run in runs.runs if not is_baseline_run(run)
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

        self.__vmax_widget = Plot()
        self.layout().replaceWidget(
            self.__ui.vmaxPlotWidget,
            self.__vmax_widget
        )

        self.__beta_widget = Plot()
        self.layout().replaceWidget(
            self.__ui.betaWidget,
            self.__beta_widget
        )

        self.__runs = runs
        self.__velocity_vs_conc = as_vel_vs_conc_series(runs)

        self.__render_plots()

    def __render_vmax_vs_km(self):
        vmax_vs_km = as_vmax_vs_km(self.__velocity_vs_conc)
        self.__vmax_widget.set_series(
            SeriesSet(
                series=[vmax_vs_km],
                x_label="1/[S]",
                y_label="1/v"
            )
        )

    def __render_plots(self):

        self.__render_vmax_vs_km()
        #exp = as_vel_vs_conc_series(self.__runs)
        #self.__plot_widget.set_series(exp)