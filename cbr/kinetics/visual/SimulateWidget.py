from PyQt5.QtCore import pyqtSlot
from typing import List, Optional

from ...core.Qt.QtWidgets import show_error
from ..compute import ComputeHandler
from ..data import *
from ..kinetics import as_conc_vs_time_series
from .FitWidgetBase import FitWidgetBase
from .Plot import PlotMeta, SeriesSet
from .PlotWidget import PlotWidget
from .Ui_SimulateWidget import Ui_SimulateWidget

class SimulateWidget(FitWidgetBase):

    def __init__(
            self,
            compute: ComputeHandler
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

        self.__fit__simulation__(
            model,
            runs.global_attributes.measurement_interval,
            runs.periods(),
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

        self.__simulate_model__(
            model,
            runs.global_attributes.measurement_interval,
            periods=runs.periods(),
            initial_concentrations=runs.concentrations()
        )

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