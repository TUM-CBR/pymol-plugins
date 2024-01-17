from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QDialog

from cbr.kinetics.kinetics import as_vel_vs_conc_series

from ..data import *
from ...core.Qt.visual.NamedTupleEditor import namedtuple_eidtor
from .LinearFitWidget import LinearFitMeta, LinearFitWidget
from .Plot import Point2d, PlotMeta
from .Ui_ParametersWizard import Ui_ParametersWizard

K_TAB_1 = "Vmax and Km"
K_TAB_2 = "Beta and Ki"

class KineticsParametersWizard(QDialog):

    def __init__(
        self
    ) -> None:
        super().__init__()

        self.__combined_velocity_vs_conc = None
        self.__ui = Ui_ParametersWizard()
        self.__ui.setupUi(self)

        self.__fit_parameters_model = namedtuple_eidtor(
            self.__ui.parametersTable,
            SubstrateInhibitionModel()
        )

        self.__vmax_widget = LinearFitWidget(
            LinearFitMeta(
                x_axis_name="1/[S]",
                y_axis_name="1/v"
            )
        )
        self.__vmax_widget.lse_model_changed.connect(self.__on_vmax_changed)
        self.__ui.initialValuesTab.insertTab(
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

        self.__ui.initialValuesTab.insertTab(
            1,
            self.__beta_widget,
            K_TAB_2
        )

        self.__ui.closeButtons.accepted.connect(self.__on_accepted)
        self.__ui.closeButtons.rejected.connect(self.__on_rejected)

    def on_runs_updated(self, runs: KineticsRuns):
        self.__combined_velocity_vs_conc = as_vel_vs_conc_series(runs, 2)
        self.__render_vmax_vs_km()

    @pyqtSlot()
    def __on_accepted(self):
        self.accept()

    @pyqtSlot()
    def __on_rejected(self):
        self.reject()

    def set_fit_parameters(self, params: SubstrateInhibitionModel):
        self.__fit_parameters_model[0] = params

    def fit_parameters(self) -> SubstrateInhibitionModel:
        params = self.__fit_parameters_model[0]
        assert params is not None, "Params should have default values"
        return params

    @pyqtSlot()
    def __on_beta_changed(self):
        model = self.__beta_widget.lse_model()
        params = self.fit_parameters()

        if model is None:
            params = params._replace(ksi = 0, beta = 0)
        else:
            beta = model.b / (1 + model.b)
            ksi = model.slope - model.slope*beta
            params = params._replace(
                beta = beta,
                ksi = ksi
            )
        self.set_fit_parameters(params)

    @pyqtSlot()
    def __on_vmax_changed(self):

        model = self.__vmax_widget.lse_model()
        params = self.fit_parameters()
        
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

        self.set_fit_parameters(params)
        self.__update_beta_and_ksi()

    def __update_beta_and_ksi(self):
        params = self.fit_parameters()
        assert params, "Params should have default values"

        velocity_vs_conc = self.__combined_velocity_vs_conc
        points = [] if velocity_vs_conc is None else velocity_vs_conc.values

        values = [
            Point2d(
                x = 1/point.x,
                y = point.y/(params.v_max - point.y)
            )
            for point in points
        ]

        self.__beta_widget.set_series(
            Series(
                PlotMeta(name=K_TAB_2),
                values
            )
        )

    def __render_vmax_vs_km(self):

        velocity_vs_conc = self.__combined_velocity_vs_conc
        points = [] if velocity_vs_conc is None else velocity_vs_conc.values

        values = [
            Point2d(
                x = 1/point.x,
                y = 1/point.y
            )
            for point in points
        ]

        self.__vmax_widget.set_series(
            Series(
                PlotMeta(name=K_TAB_1),
                values
            )
        )

        if len(values) > 0:
            self.__vmax_widget.set_x_min(values[-1].x * 0.75)
