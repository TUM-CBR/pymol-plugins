from typing import Any
from .core import *

ForEachAction = Callable[[TValue], Any]

class ForEach(ObserverBase[TValue]):

    def __init__(
        self,
        action: ForEachAction[TValue]
    ) -> None:
        super().__init__()
        self.__action = action

    def on_completed(self) -> None:
        pass

    def on_error(self, exn: Exception) -> None:
        pass

    def on_next(self, next: TValue) -> None:
        self.__action(next)