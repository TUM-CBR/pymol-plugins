from PyQt5.QtCore import pyqtSignal, pyqtSlot, QObject
from typing import List, Optional

from ..core.atomic import AtomicCounter
from ..extra.main import CbrExtraProcess, run_cbr_tools_interactive
from .data import EvalModelMetadata, Point2d, Series, SubstrateInhibitionModel

class ComputeHandler(QObject):

    MESSAGE_ID_COUNTER = AtomicCounter()
    on_model_eval_signal = pyqtSignal(object)
    on_error = pyqtSignal(object)

    def __init__(
        self,
        process : CbrExtraProcess
    ) -> None:
        super().__init__()

        self.__process = process
        self.__process.message_signal.connect(self.__on_process_message)

    def __on_eval_result(self, result: dict):
        try:
            series = Series(
                metadata = EvalModelMetadata(model_name = "partial_inhibition"),
                values = [
                    Point2d(
                        item['x'],
                        item['y']
                    )
                    for item in result['results']
                ]
            )
            self.on_model_eval_signal.emit(series)
        except IndexError:
            self.on_error.emit(ValueError("Message is not valid results"))

    def __on_value(self, value: dict):
        payload = value.get('payload')

        if payload is None:
            self.on_error.emit(ValueError("Message has no payload"))
            return

        eval_result = payload.get('eval_result')
        if eval_result is not None:
            self.__on_eval_result(eval_result)
            return

    @pyqtSlot(object)
    def __on_process_message(self, msg : dict):
        value = msg.get('value')

        if value is not None:
            self.__on_value(value)
            return

        error = msg.get('error')
        if error is not None:
            self.on_error.emit(Exception(error))
            return

    def __to_eval_message(
        self,
        model : SubstrateInhibitionModel,
        data : List[float]
    ) -> dict:
        return {
            "eval_model": {
                "model": {
                    "model_name": "partial_inhibition",
                    "model_parameters": {
                        "ksi": model.ksi,
                        "km": model.km,
                        "vmax": model.v_max,
                        "beta": model.beta
                    }
                },
                "data": data
            }
        }

    def __to_input(
        self,
        data: Optional[dict] = None,
        entity_type: str = "message"
    ):
        return {
            "uid": self.MESSAGE_ID_COUNTER.increment(),
            "entity_type": entity_type,
            "payload": data or {}
        }

    def request_model_eval(
        self,
        model : SubstrateInhibitionModel,
        data : List[float]
    ):
        message = self.__to_input(self.__to_eval_message(model, data))
        self.__process.write_json_dict(message)

def init_compute():
    return run_cbr_tools_interactive(["kinetics", "interactive"])