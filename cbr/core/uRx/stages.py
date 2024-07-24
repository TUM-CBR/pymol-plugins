from typing import Any, Generic, Optional

from .core import *

ForEachAction = Callable[[TValue], Any]
CompleteAction = Callable[[], Any]
ErrorAction = Callable[[Exception], Any]
Predicate = Callable[[TValue], bool]

class FromValues(ObservableBase[TValue]):
    def __init__(self, values: Iterable[TValue]) -> None:
        super().__init__()
        self.__values = list(values)

    def subscribe(self, os: Observer[TValue]) -> Subscription:
        for v in self.__values:
            os.on_next(v)
        os.on_completed()
        return SubscriptionEmpty()

class Filter(ObservableBase[TValue]):

    def __init__(
        self,
        os: Observable[TValue],
        predicate: Predicate[TValue],
    ) -> None:
        super().__init__()
        self.__os = os
        self.__predicate = predicate

    def subscribe(self, os: ObserverBase[TValue]) -> SubscriptionBase:

        def on_next(value: TValue):
            try:
                if self.__predicate(value):
                    os.on_next(value)
            except Exception as e:
                os.on_error(e)

        return self.__os.subscribe(
            ForEach(
                on_next,
                on_complete=lambda: os.on_completed(),
                on_error=lambda e: os.on_error(e)
            )
        )

class ExceptionWithBacktrace(Exception):
    def __init__(self, exn: Exception, backtrace: str) -> None:
        super().__init__(str(exn))
        self.__exn = exn
        self.__backtrace = backtrace

    @property
    def exn(self) -> Exception:
        return self.__exn

    @property
    def backtrace(self) -> str:
        return self.__backtrace

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
            import traceback
            self.on_error(ExceptionWithBacktrace(exn, traceback.format_exc()))

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

class Concat(ObservableBase[TValue]):

    def __init__(
        self,
        observables: Iterable[Observable[TValue]]
    ) -> None:
        super().__init__()
        self.__observables = list(observables)

    def subscribe(self, os: Observer[TValue]) -> Subscription:

        queue = list(self.__observables)
        subscription: Optional[Subscription] = None

        def on_complete():
            nonlocal queue, subscription
            
            if subscription is not None:
                subscription.dispose()
                subscription = None

            if len(queue) == 0:
                os.on_completed()
            else:
                subscription = queue.pop(0).subscribe(
                    ForEach(
                        lambda e: os.on_next(e),
                        on_error=lambda er: os.on_error(er),
                        on_complete=on_complete
                    )
                )

        def unsubscribe():
            nonlocal subscription, queue
            if subscription is not None:
                subscription.dispose()
                subscription = None
            queue = []

        on_complete()
        return SubscriptionCallable(unsubscribe)

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

TState = TypeVar('TState')

class Scan(ObservableBase[TState]):

    def __init__(
        self,
        os: Observable[TValue],
        seed: TState,
        accumulator: Callable[[TState, TValue], TState]
    ) -> None:
        super().__init__()
        self.__os = os
        self.__seed = seed
        self.__accumulator = accumulator

    def subscribe(self, os: Observer[TState]) -> Subscription:
        state = self.__seed

        def on_next(value: Any):
            nonlocal state
            state = self.__accumulator(state, value)
            os.on_next(state)

        return self.__os.subscribe(
            ForEach(
                on_next,
                on_error=os.on_error,
                on_complete=os.on_completed
            )
        )
