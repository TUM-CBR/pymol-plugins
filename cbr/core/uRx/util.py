from typing import Optional

from .core import Subscription, SubscriptionBase

class UpdatableSubscription(SubscriptionBase):

    def __init__(
        self,
        current: Optional[Subscription] = None
    ) -> None:
        super().__init__()
        self.__current = current

    def update(self, new_subscription: Subscription) -> None:

        if self.__current is not None:
            self.__current.dispose()

        self.__current = new_subscription

    def dispose(self) -> None:
        if self.__current is not None:
            self.__current.dispose()