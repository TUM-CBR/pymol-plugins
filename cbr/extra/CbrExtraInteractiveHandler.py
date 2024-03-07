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

class CbrExtraInteractiveHandlerBase(QObject):
    pass

class CbrExtraInteractiveHandler(CbrExtraInteractiveHandlerBase, Generic[TMessageIn, TMessageOut]):

    message_signal = pyqtSignal(object)

class CbrExtraInteractiveManager(QObject):

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