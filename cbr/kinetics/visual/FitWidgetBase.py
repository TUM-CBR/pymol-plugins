from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget
from typing import List

from ..compute import ComputeHandler
from ..data import EvalModelMetadata, KineticsRuns, Point2d, Series, SubstrateInhibitionModel

class FitWidgetBase(QWidget):

    def __init__(self, compute: ComputeHandler) -> None:
        super().__init__()
        self.__compute = compute
        self.__compute.on_busy_changed.connect(self.__on_busy_changed)
        self.__compute.on_model_eval_signal.connect(self.__on_model_eval)

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

    def on_model_updated(self, model: SubstrateInhibitionModel):
        self.__on_model_updated__(model)

    @pyqtSlot(object)
    def __on_model_eval(self, result: 'Series[EvalModelMetadata]'):
        self.__on_model_eval__(result)

    def on_runs_updated(self, runs: KineticsRuns):
        self.__on_runs_updated__(runs)

    def __eval_model__(
            self,
            model : SubstrateInhibitionModel,
            data: List[Point2d]
        ):
        self.__compute.request_model_eval(
            model,
            [
                point.x
                for point in data
            ]
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