from PyQt5.QtCore import pyqtSlot
from typing import List, Optional

from ...core.Qt.QtWidgets import show_error
from ...core.Qt.visual.NamedTupleEditor import namedtuple_eidtor
from ..data import *
from ..kinetics import as_conc_vs_time_series
from .FitWidgetBase import FitWidgetArgsBase, FitWidgetBase
from .Plot import PlotMeta, SeriesSet
from .PlotWidget import PlotWidget
from .Ui_SimulateWidget import Ui_SimulateWidget

class SimulateArgs(NamedTuple):
    iterations: int = 100
    steps_per_second : int = 10

class SimulateWidget(FitWidgetBase):

    def __init__(
            self,
            compute: FitWidgetArgsBase
        ) -> None:
        super().__init__(compute)

        self.__runs = None
        self.__model = None
        self.__ui = Ui_SimulateWidget()
        self.__ui.setupUi(self)

        self.__plot = PlotWidget()
        layout = self.layout()
        assert layout is not None, "This widget needs a layout."

        layout.replaceWidget(
            self.__ui.plotWidget,
            self.__plot
        )

        self.__ui.fitButton.clicked.connect(self.__on_fit_simulation)
        self.__simulate_args_model = namedtuple_eidtor(
            self.__ui.parametersWidget,
            SimulateArgs()
        )
        self.__simulate_args_model.dataChanged.connect(self.__on_parameters_changed)

    def __simulate_args(self) -> SimulateArgs:
        args = self.__simulate_args_model.current_values[0]
        assert args is not None, "The args config must have one item"
        return args

    @pyqtSlot()
    def __on_fit_simulation(self):

        runs = self.__runs
        model = self.__model

        if runs is None or model is None:
            show_error(
                self,
                "Error",
                "Data must selected before using this feature"
            )
            return
        
        args = self.__simulate_args()

        self.__fit__simulation__(
            model,
            runs.global_attributes.measurement_interval,
            runs.periods(),
            args.steps_per_second,
            args.iterations,
            dict(
                (
                    series.metadata.concentration,
                    [point.y for point in series.values]
                )
                for series in  as_conc_vs_time_series(runs)
            )
        )

    def __update_busy_state__(self, is_busy: bool):
        self.__ui.fitButton.setEnabled(not is_busy)

    def __on_model_updated__(self, model: SubstrateInhibitionModel):
        self.__model = model
        self.__update()

    def __on_runs_updated__(self, runs: KineticsRuns):
        self.__runs = runs
        scale = 1 / runs.global_attributes.concentration_units
        self.__plot.set_scaling(scale, scale)
        self.__update()

    def __update(
        self,
        model_runs: 'Optional[List[Series[SimulateModelMetadata]]]' = None
    ):

        runs = self.__runs
        if runs is None:
            return

        self.__update_plot(
            as_conc_vs_time_series(runs),
            model_runs
        )

        model = self.__model

        if model is None or model_runs is not None:
            return
        
        args = self.__simulate_args()

        self.__simulate_model__(
            model,
            runs.global_attributes.measurement_interval,
            periods = runs.periods(),
            steps_per_second = args.steps_per_second,
            initial_concentrations=runs.concentrations()
        )

    def __on_parameters_changed(self):
        self.__update()

    def __on_model_eval_simulation__(self, result: List[Series[SimulateModelMetadata]]):
        self.__update(result)        

    def __update_plot(
        self,
        experimental: List['Series[RunMetadata]'],
        model: 'Optional[List[Series[SimulateModelMetadata]]]' = None
    ):
        
        experimental_series = [
            series.update_meta(
                PlotMeta(f"Exp {series.metadata.concentration}")
            )
            for series in experimental
        ]

        if model is None:
            model_series = []
        else:
            model_series = [
                series.update_meta(
                    PlotMeta(f"Sim {series.metadata.initial_concentration}")
                )
                for series in model
            ]

        self.__plot.set_series(
            SeriesSet(
                experimental_series + model_series,
                x_label="Time",
                y_label="Concentration"
            )
        )