from typing import Any, Callable, cast, List, Optional
from PyQt5.QtCore import QObject, pyqtBoundSignal, pyqtSlot

from cbr.core.uRx.core import Observer, Subscription

from .core import *
from .dsl import *

class QtObservableBase(QObject, ObservableBase[TValue]):
            
    def __init__(self, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self.__observers: List[Observer[TValue]] = []
        self.destroyed.connect(self.__on_destroyed)

    def __with_observers__(self, fn: Callable[[Observer[TValue]], Any]) -> None:
        for observer in self.__observers:
            fn(observer)

    def __unsubscribe(self, os: Observer[TValue]) -> None:
        ix = self.__observers.index(os)

        if ix >= 0:
            self.__observers.pop(ix)

    def __on_subscribe__(self, os: ObserverBase[TValue]) -> None:
        pass

    @pyqtSlot()
    def __on_destroyed(self):
        self.__with_observers__(lambda os: os.on_completed())

    def subscribe(self, os: ObserverBase[TValue]) -> Subscription:

        self.__on_subscribe__(os)
        self.__observers.append(os)
        return SubscriptionCallable(
            lambda: self.__unsubscribe(os)
        )

class QtSignalObservable(QtObservableBase[TValue]):

    def __init__(
        self,
        parent: QObject,
        signal: pyqtBoundSignal
    ) -> None:
        super().__init__(parent)

        # With a proper programming language, we don't need two variables
        # as Some(None) != None. But Python is not such a language.
        self.__latest: Optional[TValue] = None
        signal.connect(self.__on_signal_value)
        self.__has_emission = False

    @pyqtSlot()
    def __on_signal_value(self, value: object):

        self.__has_emission = True

        # We assume that the user has verified that the type
        # of the signal is correct. They can extend the observer
        # with additional stages if more validation is needed
        self.__latest = latest = cast(TValue, value)
        self.__with_observers__(lambda os: os.on_next(latest))

    def __on_subscribe__(self, os: Observer[TValue]) -> None:

        if self.__has_emission:
            # It is not even possible to use assert to ensure consistency
            os.on_next(cast(TValue, self.__latest))

def observer_from_signal(parent: QObject, signal: pyqtBoundSignal) -> Observable[Any]:
    return QtSignalObservable(parent, signal)

def dsl_from_signal(parent: QObject, signal: pyqtBoundSignal) -> Dsl[Any]:
    return Dsl(QtSignalObservable(parent, signal))