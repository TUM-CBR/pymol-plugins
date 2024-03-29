from PyQt5.QtCore import pyqtSlot, QModelIndex
from PyQt5.QtWidgets import QDialog, QWidget
from typing import Tuple

from ...core.atomic import AtomicCounter
from ...core.Context import Context
from ...core.Qt.visual.NamedTupleEditor import namedtuple_eidtor
from ...core.Qt.QtWidgets import show_error, show_exception
from ..compute import ComputeHandler, init_compute
from ..data import *
from ..kinetics import *
from .FitWidgetBase import FitWidgetArgsBase, FitWidgetBase
from .KineticsParamtersWizard import KineticsParametersWizard
from .KineticsInput import KineticsInput
from .Plot import PlotMeta, Point2d, Series
from .SimulateWidget import SimulateWidget
from .VelocityFitWidget import VelocityFitWidget
from .Ui_KineticsOptimizer import Ui_KineticsOptimizer
from .visual import FitWidgetBase

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
        context: Context,
    ) -> None:
        super().__init__()

        self.__ui = Ui_KineticsOptimizer()
        self.__ui.setupUi(self)

        self.__runs = None
        self.__compute = ComputeHandler(init_compute())
        self.__compute.on_error.connect(self.__on_compute_error)
        self.__compute.on_model_fit_singal.connect(self.__on_fit_signal)

        self.__set_busy(False)
        self.__compute.on_busy_changed.connect(self.__on_busy_changed)

        # Create the widget for data input
        self.__data_input_widget = KineticsInput()
        data_container_layout = self.__ui.dataContainer.layout()
        assert data_container_layout is not None, "UI file does not provide a layout for dataContainer"
        data_container_layout.replaceWidget(
            self.__ui.dataWidget,
            self.__data_input_widget
        )
        self.__data_input_widget.runs_selected.connect(self.__on_runs_updated)
        initial_min = SubstrateInhibitionModel()
        initial_max = SubstrateInhibitionModel(1,1,1,1)
        args_base = FitWidgetArgsBase(self.__compute, (initial_min, initial_max))

        # Create the widget to fit by rate
        self.__velocity_widget = VelocityFitWidget(args_base)
        velocity_layout = self.__ui.fitContainer.layout()
        assert velocity_layout is not None, "UI file does not provide a layout for the velocity widget"
        velocity_layout.replaceWidget(
            self.__ui.fitWidget,
            self.__velocity_widget
        )

        # Create the widget to fit by simulation
        self.__simulation_widget = SimulateWidget(args_base)
        simulation_layout = self.__ui.simulateContainer.layout()
        assert simulation_layout is not None, "UI file does not provide a layout for the simulation widget"
        simulation_layout.replaceWidget(
            self.__ui.simulateWidget,
            self.__simulation_widget
        )

        self.__parameters_wizard = KineticsParametersWizard()
        self.__fit_parameters_model = namedtuple_eidtor(
            self.__ui.parametersTable,
            self.__parameters_wizard.fit_parameters()
        )
        self.__fit_parameters_model.dataChanged.connect(self.__on_fit_parameters_changed)
        self.__ui.wizardButton.clicked.connect(self.__on_wizard_clicked)

        self.__fit_parameters_range = namedtuple_eidtor(
            self.__ui.parametersRangeTable,
            initial_min,
            initial_max
        )
        self.__fit_parameters_range.dataChanged.connect(self.__on_ranges_changed)

    def __get_fit_widgets(self) -> Iterable[FitWidgetBase]:
        return [
            self.__velocity_widget,
            self.__simulation_widget
        ]
    
    def __get_fit_parameters_range(self):
        current_values = self.__fit_parameters_range.current_values
        min_range = current_values[0]
        max_range = current_values[1]

        assert min_range is not None and max_range is not None, "Ranges widget must have two items"
        return (
            self.__scale_model(min_range, False),
            self.__scale_model(max_range, False)
        )

    @pyqtSlot()
    def __on_ranges_changed(self):

        for widget in self.__get_fit_widgets():
            widget.on_parameters_range_changed(self.__get_fit_parameters_range())

    @pyqtSlot()
    def __on_runs_updated(self):
        self.__runs = runs = self.__data_input_widget.get_runs()

        self.__parameters_wizard.on_runs_updated(runs)

        for widget in self.__get_fit_widgets():
            widget.on_runs_updated(runs)

        self.__set_fit_parameters(self.__parameters_wizard.fit_parameters())

        # Changing the runs can potentially change the scaling, so we call
        # the "on_ranges_changed" method to indicate it.
        self.__on_ranges_changed()

    @pyqtSlot()
    def __on_wizard_clicked(self):

        runs = self.__runs

        if runs is None:
            show_error(
                self,
                "Input Error",
                "Please input the run data before using this function."
            )
            return

        self.__parameters_wizard.set_fit_parameters(self.__fit_parameters())

        if self.__parameters_wizard.exec() == QDialog.Rejected:
            return

        self.__set_fit_parameters(self.__parameters_wizard.fit_parameters())

    def __set_busy(self, busy: bool):
        self.__ui.fitProgressBar.setVisible(busy)

    @pyqtSlot()
    def __on_busy_changed(self):
        self.__set_busy(self.__compute.is_busy())

    @pyqtSlot(object)
    def __on_compute_error(self, error: Exception):
        show_exception(self, error)

    def __scale_model(self, model: SubstrateInhibitionModel, inverted: bool):
        runs = self.__runs
        factor = 1 if runs is None else runs.global_attributes.concentration_units
        if inverted:
            factor = 1 / factor
        return model.scale(factor)

    def __fit_parameters(self) -> SubstrateInhibitionModel:
        params = self.__fit_parameters_model[0]
        assert params is not None, "Params should have default values"
        return self.__scale_model(params, False)

    def __set_fit_parameters(self, params: SubstrateInhibitionModel):
        self.__fit_parameters_model[0] = self.__scale_model(params, True)

    @pyqtSlot(object)
    def __on_fit_signal(self, result: FitModelResult):
        self.__set_fit_parameters(result.model)

    @pyqtSlot(QModelIndex, QModelIndex)
    def __on_fit_parameters_changed(self, start: QModelIndex, end: QModelIndex):
        model = self.__fit_parameters()

        assert model is not None, "The model editor should have at least one model"
        self.__velocity_widget.on_model_updated(model)
        self.__simulation_widget.on_model_updated(model)