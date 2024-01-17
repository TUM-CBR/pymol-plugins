from PyQt5.QtCore import pyqtSlot
from typing import List, NamedTuple, Optional

from ...core.Qt.QtWidgets import show_error
from ...core.Qt.visual.NamedTupleEditor import namedtuple_eidtor
from ..data import *
from ..kinetics import as_vel_vs_conc_series

from .FitWidgetBase import FitWidgetArgsBase, FitWidgetBase
from .Plot import PlotMeta, SeriesSet
from .PlotWidget import PlotWidget
from .Ui_VelocityFitWidget import Ui_VelocityFitWidget

class VelocityConfiguration(NamedTuple):
    samples_to_consider: int = 2
    iterations: int = 1000

class VelocityFitWidget(FitWidgetBase):

    def __init__(self, compute: FitWidgetArgsBase) -> None:
        super().__init__(compute)
        self.__ui = Ui_VelocityFitWidget()
        self.__ui.setupUi(self)
        self.__runs = None
        self.__velocity_vs_conc = None
        self.__model = None

        self.__plot_widget = PlotWidget()
        layout = self.layout()

        assert layout is not None, "The layout should not be none. Check the ui file!"
        layout.replaceWidget(
            self.__ui.velocityPlot,
            self.__plot_widget
        )

        self.__fit_parameters_model = namedtuple_eidtor(
            self.__ui.configurationWidget,
            VelocityConfiguration()
        )

        self.__fit_parameters_model.dataChanged.connect(self.__on_parameters_changed)

        self.__ui.fitButton.clicked.connect(self.__fit_model)

    def __update_busy_state__(self, is_busy: bool):
        self.__ui.fitButton.setEnabled(not is_busy)

    def __fit_parameters(self) -> VelocityConfiguration:
        samples = self.__fit_parameters_model.current_values[0]
        assert samples is not None, "The samples configuration widget must have one element"
        return samples

    def __as_vel_vs_conc_series(self, runs: KineticsRuns) ->  'Series[RunVelocityMetadata]':
        samples = self.__fit_parameters()
        return as_vel_vs_conc_series(runs, samples.samples_to_consider)

    def __on_runs_updated__(self, runs: KineticsRuns):
        self.__runs = runs
        scale = 1 / runs.global_attributes.concentration_units
        self.__plot_widget.set_scaling(
            scale,
            scale
        )
        self.__on_runs_updated()

    def __on_parameters_changed(self):
        self.__on_runs_updated()

    def __on_runs_updated(self):

        runs = self.__runs
        if runs is None:
            self.__velocity_vs_conc = None
        else:
            self.__velocity_vs_conc = self.__as_vel_vs_conc_series(runs)
        self.__update()

    def __on_model_updated__(self, model: SubstrateInhibitionModel):
        self.__model = model
        self.__update()

    def __on_model_eval__(self, result: 'Series[EvalModelMetadata]'):
        velocity_vs_conc = self.__velocity_vs_conc

        if velocity_vs_conc is None:
            experimental = []
            model = []
        else:
            experimental = velocity_vs_conc.values
            model = result.values

        self.__update_plot(
            experimental,
            model
        )

    def __update(self):
        experimental = [] if self.__velocity_vs_conc is None \
            else self.__velocity_vs_conc.values

        self.__update_plot(experimental)
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
        
        args = self.__fit_parameters()
        self.__fit_model__(
            model,
            args.iterations,
            velocity_vs_conc.values,
        )

    def __eval_model(self):
        model = self.__model
        velocity_vs_conc = self.__velocity_vs_conc

        if model is None or velocity_vs_conc is None:
            return

        self.__eval_model__(
            model,
            velocity_vs_conc.values
        )