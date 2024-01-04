from PyQt5.QtCore import pyqtSlot, QModelIndex
from PyQt5.QtWidgets import QWidget
from typing import Any, Optional, Tuple

from ...core.atomic import AtomicCounter
from ...core.Qt.visual.NamedTupleEditor import namedtuple_eidtor
from ...core.Qt.QtWidgets import show_exception
from ..compute import ComputeHandler
from ..data import *
from ...extra.main import CbrExtraProcess
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

    MESSAGE_COUNTER = AtomicCounter()

    def __init__(
        self,
        runs: KineticsRuns,
        cbr_extra : CbrExtraProcess
    ) -> None:
        super().__init__()

        self.__ui = Ui_KineticsOptimizer()
        self.__ui.setupUi(self)

        self.__compute = ComputeHandler(cbr_extra)
        self.__compute.on_error.connect(self.__on_compute_error)
        self.__compute.on_model_eval_signal.connect(self.__on_compute_result)

        self.__fit_parameters_model = namedtuple_eidtor(
            self.__ui.parametersTable,
            SubstrateInhibitionModel()
        )
        self.__fit_parameters_model.dataChanged.connect(self.__on_fit_parameters_changed)

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
        self.__vmax_widget.lse_model_changed.connect(self.__on_vmax_changed)
        self.__ui.initalValuesWidget.insertTab(
            0,
            self.__vmax_widget,
            K_TAB_1
        )

        self.__beta_widget = LinearFitWidget(
            LinearFitMeta(
                x_axis_name="1/[S]",
                y_axis_name="v/(Vmax - v)"
            )
        )
        self.__beta_widget.lse_model_changed.connect(self.__on_beta_changed)

        self.__ui.initalValuesWidget.insertTab(
            1,
            self.__beta_widget,
            K_TAB_2
        )

        self.__runs = runs
        self.__velocity_vs_conc = as_vel_vs_conc_series(runs)

        self.__render_plots()

    @pyqtSlot(object)
    def __on_compute_error(self, error: Exception):
        show_exception(self, error)

    def __fit_parameters(self) -> SubstrateInhibitionModel:
        params = self.__fit_parameters_model[0]
        assert params is not None, "Params should have default values"
        return params

    def __set_fit_parameters(self, params: SubstrateInhibitionModel):
        self.__fit_parameters_model[0] = params

    @pyqtSlot(QModelIndex, QModelIndex)
    def __on_fit_parameters_changed(self, start, end):
        self.__render_model(model_series=None)
        self.__compute.request_model_eval(
            self.__fit_parameters(),
            [
                point.x
                for point in self.__combined_velocity_vs_conc()
            ]
        )

    @pyqtSlot(object)
    def __on_compute_result(self, result: 'Series[EvalModelMetadata]'):
        self.__render_model(model_series=result)

    @pyqtSlot()
    def __on_beta_changed(self):
        model = self.__beta_widget.lse_model()
        params = self.__fit_parameters()

        if model is None:
            params = params._replace(ksi = 0, beta = 0)
        else:
            beta = model.b / (1 + model.b)
            ksi = model.slope - model.slope*beta
            params = params._replace(
                beta = beta,
                ksi = ksi
            )
        self.__set_fit_parameters(params)

    @pyqtSlot()
    def __on_vmax_changed(self):

        model = self.__vmax_widget.lse_model()
        params = self.__fit_parameters()
        
        if model is None:
            params = params._replace(
                v_max = 0,
                km = 0
            )
        else:
            params = params._replace(
                v_max = 1/model.b,
                km = model.slope/model.b
            )

        self.__set_fit_parameters(params)
        self.__update_beta_and_ksi()

    def __combined_velocity_vs_conc(self):
        result = [
            point
            for series in self.__velocity_vs_conc
            for point in series.values
        ]
        result.sort(key=lambda point: point.x)
        return result

    def __update_beta_and_ksi(self):
        params = self.__fit_parameters()

        assert params, "Params should have default values"

        values = [
            Point2d(
                x = 1/point.x,
                y = point.y/(params.v_max - point.y)
            )
            for point in self.__combined_velocity_vs_conc()
        ]

        self.__beta_widget.set_series(
            Series(
                PlotMeta(name=K_TAB_2),
                values
            )
        )

    def __render_vmax_vs_km(self):

        values = [
            Point2d(
                x = 1/point.x,
                y = 1/point.y
            )
            for point in self.__combined_velocity_vs_conc()
        ]

        self.__vmax_widget.set_series(
            Series(
                PlotMeta(name=K_TAB_1),
                values
            )
        )

        if len(values) > 0:
            self.__vmax_widget.set_x_min(values[-1].x * 0.75)

    def __render_model(
        self,
        model_series: 'Optional[Series[Any]]'
    ):
        exp_values = self.__combined_velocity_vs_conc()
        exp_series = Series(
            metadata = PlotMeta("Experimental"),
            values = exp_values
        )

        if model_series is not None:
            series = [
                exp_series,
                model_series.update_meta(
                    PlotMeta("Fit (TF)")
                )
            ]
        else:
            series = [exp_series]

        self.__plot_widget.set_series(
            SeriesSet(
                series = series,
                x_label = "[S]",
                y_label = "v",
            )
        )

    def __render_plots(self):

        self.__render_vmax_vs_km()
        #exp = as_vel_vs_conc_series(self.__runs)
        #self.__plot_widget.set_series(exp)