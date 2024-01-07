from PyQt5.QtCore import pyqtSlot, QModelIndex
from PyQt5.QtWidgets import QDialog, QWidget
from typing import Any, Optional, Tuple

from ...core.atomic import AtomicCounter
from ...core.Qt.visual.NamedTupleEditor import namedtuple_eidtor
from ...core.Qt.QtWidgets import show_exception
from ..compute import ComputeHandler
from ..data import *
from ...extra.main import CbrExtraProcess
from ..kinetics import *
from .KineticsParamtersWizard import KineticsParametersWizard
from .Plot import Plot, PlotMeta, Point2d, Series, SeriesSet
from .SeriesModel import SeriesModel
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

class KineticsOptimizer(QWidget):

    MESSAGE_COUNTER = AtomicCounter()

    def __init__(
        self,
        runs: KineticsRuns,
        cbr_extra : CbrExtraProcess
    ) -> None:
        super().__init__()

        self.__velocity_vs_conc = as_vel_vs_conc_series(runs)
        self.__ui = Ui_KineticsOptimizer()
        self.__ui.setupUi(self)

        self.__compute = ComputeHandler(cbr_extra)
        self.__compute.on_error.connect(self.__on_compute_error)
        self.__compute.on_model_eval_signal.connect(self.__on_compute_result)
        self.__compute.on_model_fit_singal.connect(self.__on_fit_signal)

        self.__set_busy(False)
        self.__compute.on_busy_changed.connect(self.__on_busy_changed)

        self.__ui.fitModelButton.clicked.connect(self.__fit_model)

        self.__parameters_wizard = KineticsParametersWizard(self.__combined_velocity_vs_conc())
        self.__fit_parameters_model = namedtuple_eidtor(
            self.__ui.parametersTable,
            self.__parameters_wizard.fit_parameters()
        )
        self.__fit_parameters_model.dataChanged.connect(self.__on_fit_parameters_changed)
        self.__ui.wizardButton.clicked.connect(self.__on_wizard_clicked)

        self.__plot_widget = Plot()
        plotLayout = self.__ui.plotWidgetLayout.layout()
        plotLayout.replaceWidget(
            self.__ui.plotWidget,
            self.__plot_widget
        )

        self.__ui.rawDataTable.setModel(
            SeriesModel(
                self.__velocity_vs_conc,
                lambda _,c: f"[{c.metadata.concentration}]"
            )
        )

        self.__runs = runs
        self.__eval_model()

    @pyqtSlot()
    def __on_wizard_clicked(self):

        self.__parameters_wizard.set_fit_parameters(self.__fit_parameters())

        if self.__parameters_wizard.exec() == QDialog.Rejected:
            return

        self.__set_fit_parameters(self.__parameters_wizard.fit_parameters())

    def __set_busy(self, busy: bool):
        self.__ui.fitProgressBar.setVisible(busy)
        self.__ui.fitModelButton.setEnabled(not busy)

    @pyqtSlot()
    def __fit_model(self):
        self.__compute.request_model_fit(
            self.__fit_parameters(),
            self.__combined_velocity_vs_conc()
        )

    @pyqtSlot()
    def __on_busy_changed(self):
        self.__set_busy(self.__compute.is_busy())

    @pyqtSlot(object)
    def __on_compute_error(self, error: Exception):
        show_exception(self, error)

    def __fit_parameters(self) -> SubstrateInhibitionModel:
        params = self.__fit_parameters_model[0]
        assert params is not None, "Params should have default values"
        return params

    def __set_fit_parameters(self, params: SubstrateInhibitionModel):
        self.__fit_parameters_model[0] = params

    @pyqtSlot(object)
    def __on_fit_signal(self, result: FitModelResult):
        self.__set_fit_parameters(result.model)

    @pyqtSlot(QModelIndex, QModelIndex)
    def __on_fit_parameters_changed(self, start, end):
        self.__eval_model()

    def __eval_model(self):
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

    def __combined_velocity_vs_conc(self) -> List[Point2d]:
        result = [
            point
            for series in self.__velocity_vs_conc
            for point in series.values
        ]
        result.sort(key=lambda point: point.x, reverse=True)
        return result

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
                    PlotMeta("Fit")
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