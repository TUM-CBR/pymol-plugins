from PyQt5.QtCore import pyqtSignal, pyqtSlot, QObject
from typing import List, Optional, Set

from ..core.atomic import AtomicCounter
from ..extra.main import CbrExtraProcess, run_cbr_tools_interactive
from .data import EvalModelMetadata, FitModelResult, Point2d, Series, SubstrateInhibitionModel

K_UID = "uid"
K_KSI = "ksi"
K_KM = "km"
K_VMAX = "vmax"
K_BETA = "beta"
K_VALUE = "value"
K_INPUT_UIDS = "input_uids"
K_ERROR = "error"
K_MODEL_PARAMETERS = "model_parameters"
K_MODEL_NAME = "model_name"
K_PARTIAL_INHIBITION = "partial_inhibition"

class ComputeHandler(QObject):

    MESSAGE_ID_COUNTER = AtomicCounter()
    on_model_eval_signal = pyqtSignal(object)
    on_model_fit_singal = pyqtSignal(object)
    on_error = pyqtSignal(object)
    on_busy_changed = pyqtSignal()

    def __init__(
        self,
        process : CbrExtraProcess
    ) -> None:
        super().__init__()

        self.__process = process
        self.__process.message_signal.connect(self.__on_process_message)
        self.__pending_messages : Set[int] = set()

    def is_busy(self) -> bool:
        return len(self.__pending_messages) > 0

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

    def __on_fit_result(self, result: dict):
        try:
            model = result['model']
            model_name = model[K_MODEL_NAME]

            assert model_name == K_PARTIAL_INHIBITION, f"Unknown model {model_name}"
            parameters = model[K_MODEL_PARAMETERS]
            value = FitModelResult(
                model = SubstrateInhibitionModel(
                    v_max = parameters[K_VMAX],
                    beta = parameters[K_BETA],
                    km = parameters[K_KM],
                    ksi = parameters[K_KSI]
                )
            )
            self.on_model_fit_singal.emit(value)
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

        fit_result = payload.get('fit_result')
        if fit_result is not None:
            self.__on_fit_result(fit_result)
            return

    def __track_pending_ack(self, msg: dict):

        uids = []
        value_msg = msg.get(K_VALUE)
        if value_msg is not None:
            uids = value_msg[K_INPUT_UIDS]
        
        error_msg = msg.get(K_ERROR)
        if error_msg is not None:
            uids = error_msg[K_INPUT_UIDS]

        count = len(self.__pending_messages)

        if len(uids) > 0:
            self.__pending_messages.clear()

        if len(self.__pending_messages) == 0 and count > 0:
            self.on_busy_changed.emit()

    def __track_pending_send(self, msg: dict):
        self.__pending_messages.add(msg[K_UID])

        if len(self.__pending_messages) == 1:
            self.on_busy_changed.emit()

    @pyqtSlot(object)
    def __on_process_message(self, msg : dict):

        self.__track_pending_ack(msg)
        value = msg.get(K_VALUE)
        if value is not None:
            self.__on_value(value)
            return

        error = msg.get(K_ERROR)
        if error is not None:
            self.on_error.emit(Exception(error))
            return

    def __to_model_spec(self, model: SubstrateInhibitionModel) -> dict:
        return {
            "model_name": K_PARTIAL_INHIBITION,
            "model_parameters": {
                K_KSI: model.ksi,
                K_KM: model.km,
                K_VMAX: model.v_max,
                K_BETA: model.beta
            }
        }

    def __to_fit_message(
        self,
        model : SubstrateInhibitionModel,
        data : List[Point2d]
    ):
        return {
            "fit_model": {
                "model": self.__to_model_spec(model),
                "data": [
                    {"x": point.x, "y": point.y}
                    for point in data
                ]
            }
        }

    def __to_eval_message(
        self,
        model : SubstrateInhibitionModel,
        data : List[float]
    ) -> dict:
        return {
            "eval_model": {
                "model": self.__to_model_spec(model),
                "data": data
            }
        }

    def __to_input(
        self,
        data: Optional[dict] = None,
        entity_type: str = "message"
    ):
        return {
            K_UID: self.MESSAGE_ID_COUNTER.increment(),
            "entity_type": entity_type,
            "payload": data or {}
        }

    def __write_json_message(self, message: dict):
        self.__track_pending_send(message)
        self.__process.write_json_dict(message)

    def request_model_eval(
        self,
        model : SubstrateInhibitionModel,
        data : List[float]
    ):
        message = self.__to_input(self.__to_eval_message(model, data))
        self.__write_json_message(message)

    def request_model_fit(
        self,
        model : SubstrateInhibitionModel,
        data : List[Point2d]
    ):
        message = self.__to_input(self.__to_fit_message(model, data))
        self.__write_json_message(message)

def init_compute():
    return run_cbr_tools_interactive(["kinetics", "interactive"])