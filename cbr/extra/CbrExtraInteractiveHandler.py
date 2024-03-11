from abc import abstractmethod
from enum import Enum
from typing import Any, Callable, cast, Dict, Generic, List, NamedTuple, Optional, TypeVar
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot, QProcess

from ..core.uRx.core import Observable
from ..core.uRx import dsl
from ..core.uRx import qt as qtRx
from .CbrExtraProcess import CbrExtraProcess, run_cbr_tools_interactive

TMessageIn = TypeVar('TMessageIn')

TMessageOut = TypeVar('TMessageOut')

K_INPUT_UIDS = 'input_uids'
K_ERROR_CODE = 'error_code'
K_PAYLOAD = 'payload'
K_MESSAGE = 'message'

class InteractiveError(NamedTuple):
    input_uids: List[int]
    error_code: int
    payload: Optional[str]
    message: str

    @classmethod
    def parse(cls, payload: Optional[Dict[Any, Any]]) -> 'Optional[InteractiveError]':

        if payload is None:
            return None

        return InteractiveError(
            input_uids = payload[K_INPUT_UIDS],
            error_code = payload[K_ERROR_CODE],
            payload = payload.get(K_PAYLOAD),
            message = payload[K_MESSAGE]
        )

class InteractiveValue(NamedTuple):
    input_uids: List[int]
    payload: Dict[Any, Any]

    @classmethod
    def parse(cls, payload: Optional[Dict[Any, Any]]) -> 'Optional[InteractiveValue]':

        if payload is None:
            return None

        return InteractiveValue(
            input_uids = payload[K_INPUT_UIDS],
            payload = payload[K_PAYLOAD]
        )

K_UID = 'uid'
K_VALUE = 'value'
K_ERROR = 'error'

class InteractiveOutput(NamedTuple):
    uid : int
    value: Optional[InteractiveValue] = None
    error: Optional[InteractiveError] = None

    @classmethod
    def parse(cls, payload: Dict[Any, Any]) -> 'InteractiveOutput':

        return InteractiveOutput(
            uid = payload[K_UID],
            value = InteractiveValue.parse(payload.get(K_VALUE)),
            error = InteractiveError.parse(payload.get(K_ERROR))
        )
    
    def has_sender_uid(self, uid: int) -> bool:

        uids = []

        value = self.value
        if value is not None:
            uids = value.input_uids

        error = self.error
        if error is not None:
            uids = error.input_uids

        return uid in uids

class CbrMessage(Generic[TMessageIn, TMessageOut]):

    def __init__(
        self,
        input_uids: List[int],
        payload: Optional[TMessageIn],
        error_message: Optional[str]
    ):
        self.__input_uids = input_uids
        self.__payload = payload
        self.__error_message = error_message

    @property
    def payload(self) -> Optional[TMessageIn]:
        return self.__payload
    
    @property
    def error(self) -> Optional[str]:
        return self.__error_message

class SendMessageResult(NamedTuple):
    message_uid: int

class CbrExtraInteractiveHandlerBase(QObject):

    message_signal = pyqtSignal(object)

    def __init__(
        self,
        manager: 'CbrExtraInteractiveManager'
    ) -> None:
        super().__init__(manager)
        self.__manager = manager

    @abstractmethod
    def __parse_message_base__(self, payload: Dict[Any, Any]) -> Any:
        raise NotImplementedError()
    
    @abstractmethod
    def __serialize_message_base__(self, message: Any) -> Dict[Any, Any]:
        raise NotImplementedError()
    
    def __should_handle__(self, message: InteractiveOutput) -> bool:
        return True
    
    def on_message(self, message: InteractiveOutput):

        if not self.__should_handle__(message):
            return
        
        value = message.value
        input_uids = []
        if value is not None:
            input_uids = value.input_uids
            payload = self.__parse_message_base__(value.payload)
        else:
            payload = None

        error = message.error
        if error is not None:
            input_uids = error.input_uids
            error_message = error.message
        else:
            error_message = None

        result: CbrMessage[Any, Any] = CbrMessage(input_uids, payload, error_message)
        self.message_signal.emit(result)
    
    def __send_message__(self, message: Any) -> SendMessageResult:
        return self.__manager.__send_message__(self.__serialize_message_base__(message))

class MessageQueingPolicy(Enum):
    HANDLE_ALL = 0
    HANDLE_LATEST = 1

MessageParser = Callable[[Dict[Any,Any]], TMessageIn]
MessageSerializer = Callable[[TMessageOut], Dict[Any, Any]]

class CbrExtraInteractiveHandler(CbrExtraInteractiveHandlerBase, Generic[TMessageIn, TMessageOut]):

    def __init__(
        self,
        manager: 'CbrExtraInteractiveManager',
        parser: MessageParser[TMessageIn],
        serializer: MessageSerializer[TMessageOut],
        queing_policy: MessageQueingPolicy
    ) -> None:
        super().__init__(manager)
        self.__message_observable: Observable[CbrMessage[TMessageIn, TMessageOut]] = qtRx.observer_from_signal(self, self.message_signal)
        self.__parser = parser
        self.__serializer = serializer
        self.__latest_send_uid = manager.BOTTOM_UID
        self.__queing_policy = queing_policy

    def __parse_message_base__(self, payload: Dict[Any, Any]) -> Any:
        return self.__parser(payload)
    
    def __serialize_message_base__(self, message: Any) -> Dict[Any, Any]:
        return self.__serializer(message)

    def observe_message(self) -> dsl.Dsl[CbrMessage[TMessageIn, TMessageOut]]:
        return dsl.observe(self.__message_observable)
    
    def observe_values(self) -> dsl.Dsl[TMessageIn]:

        def mapping(message: CbrMessage[TMessageIn, TMessageOut]):
            if message.error is not None:
                raise Exception(message.error)
            
            assert message.payload is not None, "Message must have a payload if there is no error"
            
            return message.payload

        return self.observe_message().map(mapping)
    
    def __should_handle__(self, message: InteractiveOutput) -> bool:

        if self.__queing_policy == MessageQueingPolicy.HANDLE_LATEST:
            return message.has_sender_uid(self.__latest_send_uid)

        return True

    def send_message(self, message: TMessageOut) -> SendMessageResult:
        result = self.__send_message__(message)
        self.__latest_send_uid = result.message_uid
        return result

K_ENTITY_TYPE = 'entity_type'
K_STOP = 'stop'

class CbrProcessExit(NamedTuple):
    exit_code: int
    exit_status: QProcess.ExitStatus

class CbrExtraInteractiveManager(QObject):

    BOTTOM_UID = -1

    __message_uid = BOTTOM_UID + 1

    def __init__(
        self,
        process: CbrExtraProcess,
        parent: Optional[QObject] = None,
    ) -> None:
        super().__init__(parent)

        self.__process = process
        process.message_signal.connect(self.__on_message)
        self.__handlers: List[CbrExtraInteractiveHandlerBase] = []
        exit_obs = qtRx.observer_from_signal(
            self,
            process.finished,
            slot_args=[int, QProcess.ExitStatus],
            signal_mapper=lambda values: CbrProcessExit(*values)
        )
        error_obs = qtRx.observer_from_signal(
            self,
            process.error_signal,
            signal_mapper=lambda e: cast(Exception, e)
        )
        self.__status_obs = dsl.observe(exit_obs) \
            .merge(
                dsl.observe(error_obs).as_error(lambda e: e)
            )
        
    def observe_status(self):
        return self.__status_obs

    @pyqtSlot(object)
    def __on_message(self, message: Dict[Any, Any]):

        payload = InteractiveOutput.parse(message)

        for handler in self.__handlers:
            handler.on_message(payload)

    def __send_message__(self, payload: Dict[Any, Any]) -> SendMessageResult:
        uid = self.__message_uid
        self.__message_uid += 1

        self.__process.write_json_dict({
            K_UID: uid,
            K_ENTITY_TYPE: K_MESSAGE,
            K_PAYLOAD: payload
        })

        return SendMessageResult(message_uid=uid)
    
    def message_handler(
        self,
        parser: MessageParser[TMessageIn],
        serializer: MessageSerializer[TMessageOut],
        queing_policy: MessageQueingPolicy = MessageQueingPolicy.HANDLE_ALL
    ) -> CbrExtraInteractiveHandler[TMessageIn, TMessageOut]:
        
        handler = CbrExtraInteractiveHandler(
            self,
            parser,
            serializer,
            queing_policy
        )

        self.__handlers.append(handler)

        return handler
    
def run_interactive(
    args : List[str],
    parent: Optional[QObject] = None
) -> CbrExtraInteractiveManager:

    return CbrExtraInteractiveManager(
        run_cbr_tools_interactive(args, parent),
        parent
    )