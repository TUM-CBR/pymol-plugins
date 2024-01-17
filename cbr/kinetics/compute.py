from PyQt5.QtCore import pyqtSignal, pyqtSlot, QObject
from typing import Any, Dict, List, Optional, Set, Tuple

from ..core.atomic import AtomicCounter
from ..extra.main import CbrExtraProcess, run_cbr_tools_interactive
from .data import *

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

#FitArgs
K_FIT_ITERATIONS = "iterations"
K_FIT_RANGE = "fit_range"

# Simulation spec keys
K_PERIODS = "periods"
K_INTERVAL = "interval"
K_STEPS_PER_SECOND = "reps_per_unit"

# Simulation Result keys
K_SIMULATION_SPEC = "simulation_spec"
K_SIMULATION_RESULT = "results"

# Interactive Output keys
K_SIMULATE_RESULT = "simulate_result"

# Fit simulation args
K_FIT_SIMULATION_MODEL = "model"
K_FIT_SIMULATION_SPEC = "simulation_spec"
K_FIT_SIMULATION_DATA = "data"
K_FIT_SIMULATION_ITERATIONS = "iterations"

# Interactive input keys
K_FIT_SIMULATION = "fit_simulation"

# Fit Range
K_RANGE_KSI = "ksi"
K_RANGE_KM = "km"
K_RANGE_VMAX = "vmax"
K_RANGE_BETA = "beta"

class ComputeHandler(QObject):

    MESSAGE_ID_COUNTER = AtomicCounter()
    on_model_eval_signal = pyqtSignal(object)
    on_model_eval_simulation_signal = pyqtSignal(object)
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
        self.__process.error_signal.connect(self.on_error)
        self.__pending_messages : Set[int] = set()

    def is_busy(self) -> bool:
        return len(self.__pending_messages) > 0

    def __on_eval_result(self, result: Dict[Any, Any]):
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

    def __on_fit_result(self, result: Dict[Any, Any]):
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

    def __on_simulate_result(self, result: Dict[Any, Any]):
        try:
            spec_dict = result[K_SIMULATION_SPEC]
            spec = SimulationSpec(
                periods=spec_dict[K_PERIODS],
                interval=spec_dict[K_INTERVAL]
            )
            results = [
                Series(
                    metadata = SimulateModelMetadata(
                        model_name="partial_inhibition",
                        initial_concentration=conc,
                        spec = spec
                    ),
                    values = [
                        Point2d(time, value)
                        for time,value in zip(spec.times(), values)
                    ]
                )
                for (conc, values) in result[K_SIMULATION_RESULT].items()
            ]
            self.on_model_eval_simulation_signal.emit(results)
        except IndexError:
            self.on_error.emit(ValueError("Message is not valid results"))

    def __on_value(self, value: Dict[Any, Any]):
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
        
        simulate_result = payload.get(K_SIMULATE_RESULT)
        if simulate_result is not None:
            self.__on_simulate_result(simulate_result)
            return

    def __track_pending_ack(self, msg: Dict[Any, Any]):

        uids : List[int] = []
        value_msg = msg.get(K_VALUE)
        if value_msg is not None:
            uids = value_msg[K_INPUT_UIDS]
        
        error_msg = msg.get(K_ERROR)
        if error_msg is not None:
            uids = error_msg[K_INPUT_UIDS]

        count = len(self.__pending_messages)

        for uid in uids:
            self.__pending_messages.remove(uid)

        if len(self.__pending_messages) == 0 and count > 0:
            self.on_busy_changed.emit()

    def __track_pending_send(self, msg: Dict[Any, Any]):
        self.__pending_messages.add(msg[K_UID])

        if len(self.__pending_messages) == 1:
            self.on_busy_changed.emit()

    @pyqtSlot(object)
    def __on_process_message(self, msg : Dict[Any, Any]):

        self.__track_pending_ack(msg)
        value = msg.get(K_VALUE)
        if value is not None:
            self.__on_value(value)
            return

        error = msg.get(K_ERROR)
        if error is not None:
            self.on_error.emit(Exception(error))
            return

    def __to_model_spec(self, model: SubstrateInhibitionModel) -> Dict[Any, Any]:
        return {
            "model_name": K_PARTIAL_INHIBITION,
            "model_parameters": {
                K_KSI: model.ksi,
                K_KM: model.km,
                K_VMAX: model.v_max,
                K_BETA: model.beta
            }
        }
    
    def __to_fit_range(
        self,
        fit_range: SubstrateInhibitionModelFitRange,    
    ) -> Dict[Any, Any]:
        min_values, max_values = fit_range
        return {
            K_RANGE_KSI: [min_values.ksi, max_values.ksi],
            K_RANGE_KM: [min_values.km, max_values.km],
            K_RANGE_VMAX: [min_values.v_max, max_values.v_max],
            K_RANGE_BETA: [min_values.beta, max_values.beta]
        }

    def __to_fit_message(
        self,
        model : SubstrateInhibitionModel,
        fit_range: SubstrateInhibitionModelFitRange,
        iterations : int,
        data : List[Point2d]
    ):
        return {
            "fit_model": {
                "model": self.__to_model_spec(model),
                "data": [
                    {"x": point.x, "y": point.y}
                    for point in data
                ],
                K_FIT_ITERATIONS: iterations,
                K_FIT_RANGE: self.__to_fit_range(fit_range)
            }
        }

    def __to_eval_message(
        self,
        model : SubstrateInhibitionModel,
        data : List[float]
    ) -> Dict[Any, Any]:
        return {
            "eval_model": {
                "model": self.__to_model_spec(model),
                "data": data
            }
        }
    
    def __to_simulation_spec(
        self,
        interval: int,
        periods : int,
        steps_per_second: int
    ):
        return {
            K_INTERVAL: interval,
            K_PERIODS: periods,
            K_STEPS_PER_SECOND: steps_per_second
        }
    
    def __to_simulate_message(
        self,
        model: SubstrateInhibitionModel,
        interval: int,
        periods: int,
        steps_per_second: int,
        initial_concentrations : List[float]
    ):
        return {
            "simulate_model": {
                "model": self.__to_model_spec(model),
                K_SIMULATION_SPEC: self.__to_simulation_spec(interval, periods, steps_per_second),
                "initial_concentrations": initial_concentrations
            }
        }
    
    def __to_fit_simulation_message(
            self,
            model: SubstrateInhibitionModel,
            fit_range: SubstrateInhibitionModelFitRange,
            interval: int,
            periods: int,
            steps_per_second: int,
            iterations: int,
            data: Dict[float, List[float]]
    ):
        return {
            K_FIT_SIMULATION: {
                K_FIT_SIMULATION_MODEL: self.__to_model_spec(model),
                K_FIT_SIMULATION_SPEC: self.__to_simulation_spec(interval, periods, steps_per_second),
                K_FIT_SIMULATION_DATA: data,
                K_FIT_SIMULATION_ITERATIONS: iterations,
                K_FIT_RANGE: self.__to_fit_range(fit_range)
            }
        }

    def __to_input(
        self,
        data: Optional[Dict[Any, Any]] = None,
        entity_type: str = "message"
    ):
        return {
            K_UID: self.MESSAGE_ID_COUNTER.increment(),
            "entity_type": entity_type,
            "payload": data or {}
        }

    def __write_json_message(self, message: Dict[Any, Any]):
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
        fit_range: SubstrateInhibitionModelFitRange,
        iterations: int,
        data : List[Point2d]
    ):
        message = self.__to_input(self.__to_fit_message(model, fit_range, iterations, data))
        self.__write_json_message(message)

    def request_model_simulate(
            self,
            model: SubstrateInhibitionModel,
            intervals: int,
            periods: int,
            steps_per_second: int,
            initial_concentrations : List[float]
    ):
        message = self.__to_input(
            self.__to_simulate_message(
                model,
                intervals,
                periods,
                steps_per_second,
                initial_concentrations
            )
        )
        self.__write_json_message(message)

    def request_fit_by_simulation(
        self,
        model: SubstrateInhibitionModel,
        fit_range: Tuple[SubstrateInhibitionModel, SubstrateInhibitionModel],
        interval: int,
        periods: int,
        steps_per_second: int,
        iterations: int,
        data : Dict[float, List[float]]
    ):
        message = self.__to_input(
            self.__to_fit_simulation_message(
                model,
                fit_range,
                interval,
                periods,
                steps_per_second,
                iterations,
                data
            )
        )
        self.__write_json_message(message)        

def init_compute():
    return run_cbr_tools_interactive(["kinetics", "interactive"])