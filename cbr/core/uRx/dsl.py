from typing import cast, Union

from cbr.core.uRx.core import Subscription

from .core import *
from .stages import *

TMappedValue = TypeVar('TMappedValue')
TMergeLeft = TypeVar('TMergeLeft')
TMergeRight = TypeVar('TMergeRight')

MergedValues = Union[TMergeLeft, TMergeRight]
MergedObservable = Observable[MergedValues[TMergeLeft, TMergeRight]]
ErrorMapper = Callable[[TValue], Exception]

class Dsl(ObservableBase[TValue], Generic[TValue]):

    def __init__(
        self,
        obs: Observable[TValue]
    ) -> None:
        self.__obs = obs

    def for_each(
        self,
        action: ForEachAction[TValue],
        on_complete: Optional[CompleteAction] = None,
        on_error: Optional[ErrorAction] = None
    ) -> Subscription:
        return self.__obs.subscribe(
            ForEach(
                action,
                on_complete,
                on_error
            )
        )
    
    def map(self, mapping: Mapping[TValue, TMappedValue]) -> 'Dsl[TMappedValue]':
        return Dsl(Map(self.__obs, mapping))
    
    def merge_union(self, right: Observable[TMergeRight]) -> 'Dsl[MergedValues[TValue, TMergeRight]]':

        # Python is a dodgey language, so Unions are not real types
        c_left : MergedObservable[TValue, TMergeRight] = cast(MergedObservable[TValue, TMergeRight], self.__obs)
        c_right : MergedObservable[TValue, TMergeRight] = cast(MergedObservable[TValue, TMergeRight], right)

        m: MergedObservable[TValue, TMergeRight] = Merge([c_left, c_right])
        return Dsl(m)
    
    def merge(self, *args: Observable[TValue]) -> 'Dsl[TValue]':
        return Dsl(Merge([self.__obs] + list(args)))
    
    def subscribe(self, os: ObserverBase[TValue]) -> SubscriptionBase:
        return self.__obs.subscribe(os)
    
    def catch(self, catch: ErrorAction):
        return Dsl(Catch(self.__obs, catch))
    
    def as_error(self, mapper: ErrorMapper[TValue]) -> 'Dsl[Any]':

        def error_mapper(v: TValue) -> Any:
            raise mapper(v)
        
        return self.map(error_mapper)
    
def observe(os: Observable[TValue]) -> Dsl[TValue]:
    return Dsl(os)