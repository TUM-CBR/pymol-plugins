from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget
from typing import Dict, List

from ..compute import ComputeHandler
from ..data import *

class FitWidgetBase(QWidget):

    def __init__(self, compute: ComputeHandler) -> None:
        super().__init__()
        self.__compute = compute
        self.__compute.on_busy_changed.connect(self.__on_busy_changed)
        self.__compute.on_model_eval_signal.connect(self.__on_model_eval)
        self.__compute.on_model_eval_simulation_signal.connect(self.__on_model_eval_simulation)

    def __compute__(self):
        return self.__compute
    
    @pyqtSlot()
    def __on_busy_changed(self):
        is_busy = self.__compute.is_busy()
        self.__update_busy_state__(is_busy)
    
    def __update_busy_state__(self, is_busy: bool):
        pass

    def __on_runs_updated__(self, runs: KineticsRuns):
        pass

    def __on_model_eval__(self, result: 'Series[EvalModelMetadata]'):
        pass

    def __on_model_updated__(self, model: SubstrateInhibitionModel):
        pass

    def __on_model_eval_simulation__(self, result: 'List[Series[SimulateModelMetadata]]'):
        pass

    def on_model_updated(self, model: SubstrateInhibitionModel):
        self.__on_model_updated__(model)

    @pyqtSlot(object)
    def __on_model_eval_simulation(self, result: 'List[Series[SimulateModelMetadata]]'):
        self.__on_model_eval_simulation__(result)

    @pyqtSlot(object)
    def __on_model_eval(self, result: 'Series[EvalModelMetadata]'):
        self.__on_model_eval__(result)

    def on_runs_updated(self, runs: KineticsRuns):
        self.__on_runs_updated__(runs)

    def __eval_model__(
            self,
            model: SubstrateInhibitionModel,
            data: List[Point2d]
        ):
        self.__compute.request_model_eval(
            model,
            [
                point.x
                for point in data
            ]
        )

    def __simulate_model__(
            self,
            model: SubstrateInhibitionModel,
            interval: int,
            periods: int,
            initial_concentrations: List[float]
    ):
        self.__compute.request_model_simulate(
            model,
            interval,
            periods,
            initial_concentrations
        )

    def __fit_model__(
            self,
            model: SubstrateInhibitionModel,
            data: List[Point2d]
    ):
        self.__compute.request_model_fit(
            model,
            data
        )

    def __fit__simulation__(
        self,
        model: SubstrateInhibitionModel,
        interval: int,
        periods: int,
        data: Dict[float, List[float]]
    ):
        self.__compute.request_fit_by_simulation(
            model,
            interval,
            periods,
            data
        )