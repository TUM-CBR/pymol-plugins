
from abc import ABC, abstractmethod
from typing import Callable, Generic, Optional, TypeVar

TValue = TypeVar('TValue')

class ObserverBase(ABC, Generic[TValue]):

    @abstractmethod
    def on_next(self, next: TValue) -> None:
        raise NotImplementedError()
    
    @abstractmethod
    def on_completed(self) -> None:
        raise NotImplementedError()
    
    @abstractmethod
    def on_error(self, exn: Exception) -> None:
        raise NotImplementedError()

Observer = ObserverBase[TValue]

class SubscriptionBase(ABC):

    @abstractmethod
    def dispose(self) -> None:
        raise NotImplementedError()

Subscription = SubscriptionBase

class SubscriptionCallable(SubscriptionBase):

    def __init__(self, unsubcribe: Callable[[], None]) -> None:
        super().__init__()
        self.__unsubscribe : Optional[Callable[[], None]] = unsubcribe

    def dispose(self) -> None:

        unsubscribe = self.__unsubscribe

        if unsubscribe is None:
            return
        
        self.__unsubscribe = None
        unsubscribe()

class ObservableBase(ABC, Generic[TValue]):

    @abstractmethod
    def subscribe(self, os: Observer[TValue]) -> Subscription:
        raise NotImplementedError()
    
Observable = ObservableBase[TValue]