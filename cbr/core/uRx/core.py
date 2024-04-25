
from abc import ABC, abstractmethod
from typing import Callable, Generic, Iterable, List, Optional, TypeVar

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

class SubscriptionComposite(SubscriptionBase):

    def __init__(
        self,
        nested: Iterable[SubscriptionBase]
    ) -> None:
        super().__init__()
        self.__subscriptions: List[SubscriptionBase] = list(nested)

    def dispose(self) -> None:
        exn = None
        for sub in self.__subscriptions:
            try:
                sub.dispose()
            except Exception as e:
                exn = e

        self.__subscriptions = []

        if exn is not None:
            raise exn
        
    def append(self, *subscriptions: SubscriptionBase) -> None:

        for subscription in subscriptions:
            self.__subscriptions.append(subscription)

class ObservableBase(Generic[TValue]):

    @abstractmethod
    def subscribe(self, os: Observer[TValue]) -> Subscription:
        raise NotImplementedError()
    
Observable = ObservableBase[TValue]