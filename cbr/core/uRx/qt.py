from typing import Any, Callable, Sequence, Type, cast, List, Optional
from PyQt5.QtCore import QObject, pyqtBoundSignal, pyqtSlot

from . import dsl
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

    @abstractmethod
    def __on_subscribe__(self, os: Observer[TValue]) -> None:
        raise NotImplementedError()

    @pyqtSlot()
    def __on_destroyed(self):
        self.__with_observers__(lambda os: os.on_completed())
        self.__observers = []

    def subscribe(self, os: ObserverBase[TValue]) -> Subscription:

        self.__on_subscribe__(os)
        self.__observers.append(os)
        return SubscriptionCallable(
            lambda: self.__unsubscribe(os)
        )
    
    def observe(self) -> Dsl[TValue]:
        return dsl.observe(self)

class QtSignalObservable(QtObservableBase[TValue]):

    def __init__(
        self,
        parent: QObject,
        hot: bool = True
    ) -> None:
        super().__init__(parent)

        # With a proper programming language, we don't need two variables
        # as Some(None) != None. But Python is not such a language.
        self.__latest: Optional[TValue] = None
        self.__has_emission = False
        self.__hot = hot

    def __on_signal_value__(self, value: object):

        self.__has_emission = True

        # We assume that the user has verified that the type
        # of the signal is correct. They can extend the observer
        # with additional stages if more validation is needed
        self.__latest = latest = cast(TValue, value)
        self.__with_observers__(lambda os: os.on_next(latest))

    def __on_subscribe__(self, os: Observer[TValue]) -> None:

        if self.__hot and self.__has_emission:
            # It is not even possible to use assert to ensure consistency
            os.on_next(cast(TValue, self.__latest))

def default_mapper(args: Sequence[Any]) -> Any:
    return next(args.__iter__())

def observer_from_signal(
    parent: QObject,
    signal: pyqtBoundSignal,
    slot_args: Optional[List[Type[Any]]] = None,
    signal_mapper: Callable[[Sequence[Any]], TValue] = default_mapper,
    hot: bool = True
) -> QtObservableBase[TValue]:
    
    if slot_args is None:
        slot_args = [object]

    class QtCustomObservable(QtSignalObservable[Any]):

        def __init__(
            self,
            parent: QObject,
            signal: pyqtBoundSignal,
            mapper: Callable[[Sequence[Any]], TValue]
        ) -> None:
            super().__init__(parent, hot=hot)
            signal.connect(self.__on_signal)
            self.__mapper = mapper

        @pyqtSlot(*slot_args)
        def __on_signal(self, *args: Any):
            self.__on_signal_value__(self.__mapper(args))

    return QtCustomObservable(parent, signal, signal_mapper)

def signal_mapper0(*_args: Any) -> tuple[()]:
    return ()

def observer_from_signal0(
    parent: QObject,
    signal: pyqtBoundSignal,
    hot: bool = True
) -> QtObservableBase[tuple[()]]:
    return observer_from_signal(
        parent,
        signal,
        slot_args=[],
        signal_mapper=signal_mapper0,
        hot=hot
    )

def dsl_from_signal(parent: QObject, signal: pyqtBoundSignal) -> Dsl[Any]:
    return Dsl(observer_from_signal(parent, signal))
