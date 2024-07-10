from typing import NamedTuple, cast, Union

from .core import *
from .stages import *

TMappedValue = TypeVar('TMappedValue')
TMergeLeft = TypeVar('TMergeLeft')
TMergeRight = TypeVar('TMergeRight')

MergedValues = Union[TMergeLeft, TMergeRight]
MergedObservable = Observable[MergedValues[TMergeLeft, TMergeRight]]
ErrorMapper = Callable[[TValue], Exception]

SubscribeCallback = Callable[[SubscriptionBase], Any]

class DslContext(NamedTuple):

    subscribe_callback: Optional[List[SubscribeCallback]] = None

    def on_subscribe(self, s: SubscriptionBase):

        if self.subscribe_callback is None:
            return
        
        for cb in self.subscribe_callback:
            cb(s)

    def add_callback(self, *cb: SubscribeCallback) -> 'DslContext':

        if self.subscribe_callback is None:
            prev = []
        else:
            prev = self.subscribe_callback

        return self._replace(subscribe_callback = prev + list(cb))

class Dsl(ObservableBase[TValue], Generic[TValue]):

    def __init__(
        self,
        obs: Observable[TValue],
        context: DslContext
    ) -> None:
        self.__obs = obs
        self.__context = context

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
        return Dsl(Map(self.__obs, mapping), self.__context)
    
    def merge_union(self, right: Observable[TMergeRight]) -> 'Dsl[MergedValues[TValue, TMergeRight]]':

        # Python is a dodgey language, so Unions are not real types
        c_left : MergedObservable[TValue, TMergeRight] = cast(MergedObservable[TValue, TMergeRight], self.__obs)
        c_right : MergedObservable[TValue, TMergeRight] = cast(MergedObservable[TValue, TMergeRight], right)

        m: MergedObservable[TValue, TMergeRight] = Merge([c_left, c_right])
        return Dsl(m, self.__context)
    
    def merge(self, *args: Observable[TValue]) -> 'Dsl[TValue]':
        return Dsl(Merge([self.__obs] + list(args)), self.__context)
    
    def with_subscribe_callback(self, *cbs: SubscribeCallback) -> 'Dsl[TValue]':
        return Dsl(
            self.__obs,
            self.__context.add_callback(*cbs)
        )

    def subscribe(self, os: ObserverBase[TValue]) -> SubscriptionBase:
        s = self.__obs.subscribe(os)
        self.__context.on_subscribe(s)

        return s
    
    def catch(self, catch: ErrorAction):
        return Dsl(Catch(self.__obs, catch), self.__context)
    
    def as_error(self, mapper: ErrorMapper[TValue]) -> 'Dsl[Any]':

        def error_mapper(v: TValue) -> Any:
            raise mapper(v)
        
        return self.map(error_mapper)
    
def observe(os: Observable[TValue]) -> Dsl[TValue]:
    return Dsl(os, DslContext())
