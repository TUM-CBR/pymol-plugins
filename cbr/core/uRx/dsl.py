from .core import *
from .stages import *

class Dsl(Generic[TValue]):

    def __init__(
        self,
        obs: Observable[TValue]
    ) -> None:
        self.__obs = obs

    def for_each(self, action: ForEachAction[TValue]) -> Subscription:
        return self.__obs.subscribe(
            ForEach(action)
        )
    
def observe(os: Observable[TValue]) -> Dsl[TValue]:
    return Dsl(os)