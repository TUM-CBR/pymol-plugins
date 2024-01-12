from typing import NamedTuple

from ..compute import ComputeHandler
from ..data import KineticsRuns, SubstrateInhibitionModel
from ..kinetics import as_vel_vs_conc_series

from .FitWidgetBase import FitWidgetBase
from .Ui_VelocityFitWidget import Ui_VelocityFitWidget

class VelocityConfiguration(NamedTuple):
    samples_to_consider: int

class VelocityFitWidget(FitWidgetBase):

    def __init__(self, compute: ComputeHandler) -> None:
        super().__init__(compute)
        self.__ui = Ui_VelocityFitWidget()
        self.__velocity_vs_conc = None
        self.__model = None

    def __update_busy_state__(self, is_busy: bool):
        self.__ui.fitButton.setEnabled(not is_busy)

    def on_runs_updated(self, runs: KineticsRuns):
        self.__velocity_vs_conc = as_vel_vs_conc_series(runs)
        self.__update()

    def on_model_updated(self, model: SubstrateInhibitionModel):
        self.__model = model
        self.__update()

    def __update(self):
        pass