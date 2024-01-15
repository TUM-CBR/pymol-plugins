from PyQt5.QtCore import pyqtSlot
from typing import List, NamedTuple, Optional

from ...core.Qt.QtWidgets import show_error
from ..compute import ComputeHandler
from ..data import EvalModelMetadata, KineticsRuns, Point2d, Series, SubstrateInhibitionModel
from ..kinetics import as_vel_vs_conc_series, combined_velocity_vs_conc

from .FitWidgetBase import FitWidgetBase
from .Plot import PlotMeta, SeriesSet
from .PlotWidget import PlotWidget
from .Ui_VelocityFitWidget import Ui_VelocityFitWidget

class VelocityConfiguration(NamedTuple):
    samples_to_consider: int

class VelocityFitWidget(FitWidgetBase):

    def __init__(self, compute: ComputeHandler) -> None:
        super().__init__(compute)
        self.__ui = Ui_VelocityFitWidget()
        self.__ui.setupUi(self)
        self.__velocity_vs_conc = None
        self.__model = None

        self.__plot_widget = PlotWidget()
        layout = self.layout()

        assert layout is not None, "The layout should not be none. Check the ui file!"
        layout.replaceWidget(
            self.__ui.velocityPlot,
            self.__plot_widget
        )

        self.__ui.fitButton.clicked.connect(self.__fit_model)

    def __update_busy_state__(self, is_busy: bool):
        self.__ui.fitButton.setEnabled(not is_busy)

    def __on_runs_updated__(self, runs: KineticsRuns):
        self.__velocity_vs_conc = as_vel_vs_conc_series(runs)
        self.__update()

    def __on_model_updated__(self, model: SubstrateInhibitionModel):
        self.__model = model
        self.__update()

    def __on_model_eval__(self, result: 'Series[EvalModelMetadata]'):
        self.__update_plot(
            combined_velocity_vs_conc(self.__velocity_vs_conc or []),
            result.values
        )

    def __update(self):
        self.__update_plot(combined_velocity_vs_conc(self.__velocity_vs_conc or []))
        self.__eval_model()

    def __update_plot(
            self,
            experimental : List[Point2d],
            model : Optional[List[Point2d]] = None
        ):

        exp_series = Series(metadata=PlotMeta("Experimental"), values = experimental)

        model_series = \
            None if model is None \
            else Series(metadata=PlotMeta("Model"), values = model)

        plot_series = [
            series
            for series in [exp_series, model_series] if series is not None
        ]

        series_set = SeriesSet(
            series = plot_series,
            x_label = "[S]",
            y_label = "v"
        )

        self.__plot_widget.set_series(
            series_set
        )

    @pyqtSlot()
    def __fit_model(self):

        velocity_vs_conc = self.__velocity_vs_conc
        model = self.__model
        if velocity_vs_conc is None or model is None:
            show_error(
                self,
                "Error",
                "You must provide data before fitting the model."
            )
            return

        self.__fit_model__(
            model,
            combined_velocity_vs_conc(velocity_vs_conc)
        )

    def __eval_model(self):
        model = self.__model

        if model is None:
            return

        self.__eval_model__(
            model,
            combined_velocity_vs_conc(
                self.__velocity_vs_conc or []
            )
        )