from abc import ABC, abstractmethod
from enum import Enum
from typing import Any, Dict, Generic, List, NamedTuple, Optional, TypeVar
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot

from .CbrExtraProcess import CbrExtraProcess

TMessageIn = TypeVar('TMessageIn')

TMessageOut = TypeVar('TMessageOut')

class InteractiveError(NamedTuple):
    input_uids: List[int]
    error_code: int
    payload: Optional[str]
    message: str

class InteractiveValue(NamedTuple):
    input_uids: List[int]
    payload: Dict[Any, Any]

class InteractiveOutput(NamedTuple):
    uid : int
    value: Optional[InteractiveValue] = None
    error: Optional[InteractiveError] = None

    @classmethod
    def parse(cls, payload: Dict[Any, Any]) -> 'InteractiveOutput':
        pass

class CbrMessage(Generic[TMessageIn, TMessageOut]):
    pass

class SendMessageResult(NamedTuple):
    message_uid: int

class CbrExtraInteractiveHandlerBase(QObject, ABC):

    def __init__(
        self,
        manager: 'CbrExtraInteractiveManager',
        parent: Optional[QObject] = None
    ) -> None:
        super().__init__(parent)

        self.__manager = manager

    @abstractmethod
    def __parse_message_base__(self, payload: Dict[Any, Any]) -> Any:
        raise NotImplementedError()
    
    @abstractmethod
    def __serialize_message_base__(self, message: Any) -> Dict[Any, Any]:
        raise NotImplementedError()
    
    def __should_handle__(self, message: InteractiveOutput) -> bool:
        return True
    
    def __send_message__(self, message: Any) -> SendMessageResult:
        return self.__manager.__send_message__(self.__serialize_message_base__(message))

class MessageQueingPolicy(Enum):
    HANDLE_ALL = 0
    HANDLE_LATEST = 1

class CbrExtraInteractiveHandler(CbrExtraInteractiveHandlerBase, Generic[TMessageIn, TMessageOut]):

    message_signal = pyqtSignal(object)

K_UID = 'uid'
K_ENTITY_TYPE = 'entity_type'
K_PAYLOAD = 'payload'
K_MESSAGE = 'message'
K_STOP = 'stop'

class CbrExtraInteractiveManager(QObject):

    __message_uid = 0

    def __init__(
        self,
        process: CbrExtraProcess,
        parent: Optional[QObject] = None,
    ) -> None:
        super().__init__(parent)

        self.__process = process
        process.message_signal.connect(self.__on_message)
        self.__handlers: List[CbrExtraInteractiveHandlerBase] = []

    @pyqtSlot()
    def __on_message(self, message: Dict[Any, Any]):
        pass

    def __send_message__(self, payload: Dict[Any, Any]) -> SendMessageResult:
        uid = self.__message_uid
        self.__message_uid += 1

        self.__process.write_json_dict({
            K_UID: uid,
            K_ENTITY_TYPE: K_MESSAGE,
            K_PAYLOAD: payload
        })

        return SendMessageResult(message_uid=uid)