from typing import Any, Generic

from .core import *

ForEachAction = Callable[[TValue], Any]
CompleteAction = Callable[[], Any]
ErrorAction = Callable[[Exception], Any]

class ForEach(ObserverBase[TValue]):

    def __init__(
        self,
        action: ForEachAction[TValue],
        on_complete: Optional[CompleteAction] = None,
        on_error: Optional[ErrorAction] = None
    ) -> None:
        super().__init__()
        self.__action = action
        self.__on_complete = on_complete
        self.__on_error = on_error

    def on_completed(self) -> None:
        if self.__on_complete is not None:
            self.__on_complete()

    def on_error(self, exn: Exception) -> None:
        if self.__on_error is not None:
            return self.__on_error(exn)
        else:
            raise exn

    def on_next(self, next: TValue) -> None:

        try:
            self.__action(next)
        except Exception as exn:
            self.on_error(exn)

TMapIn = TypeVar('TMapIn')
Mapping = Callable[[TMapIn], TValue]

class Map(ObservableBase[TValue], Generic[TMapIn, TValue]):

    def __init__(
        self,
        os: Observable[TMapIn],
        mapping: Mapping[TMapIn, TValue]
    ):
        self.__os = os
        self.__mapping = mapping

    def subscribe(self, os: ObserverBase[TValue]) -> SubscriptionBase:

        def on_next(value: TMapIn):
            try:
                os.on_next(self.__mapping(value))
            except Exception as e:
                os.on_error(e)

        return self.__os.subscribe(
            ForEach(
                on_next,
                on_complete=lambda: os.on_completed(),
                on_error=lambda e: os.on_error(e)
            )
        )

class Merge(ObservableBase[TValue]):
    def __init__(
        self,
        observables: Iterable[Observable[TValue]]
    ) -> None:
        super().__init__()
        self.__observables = list(observables)

    def subscribe(self, os: Observer[TValue]) -> Subscription:

        total_subscriptions = len(self.__observables)

        def on_complete():
            nonlocal os, total_subscriptions

            total_subscriptions -= 1

            if total_subscriptions < 0:
                raise Exception("The 'on_complete' method must only be called once.")
            elif total_subscriptions == 0:
                os.on_completed()

        subscriptions = [
            obs.subscribe(
                ForEach(
                    lambda e: os.on_next(e),
                    on_error=lambda er: os.on_error(er),
                    on_complete=on_complete
                )
            )
            for obs in self.__observables
        ]

        return SubscriptionComposite(subscriptions)

class Catch(ObservableBase[TValue]):

    def __init__(
        self,
        os: Observable[TValue],
        on_error: ErrorAction
    ) -> None:
        super().__init__()
        self.__obs = os
        self.__on_error = on_error

    def subscribe(self, os: Observer[TValue]) -> Subscription:
        return self.__obs.subscribe(
            ForEach(
                os.on_next,
                on_error=self.__on_error,
                on_complete=os.on_completed
            )
        )
    