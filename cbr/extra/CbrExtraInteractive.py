from PyQt5.QtCore import QObject
from typing import Generic, List, Optional

from ..core.uRx.core import SubscriptionBase, SubscriptionComposite
from ..core.uRx.dsl import Dsl
from .CbrExtraInteractiveHandler import CbrExtraInteractiveHandler, CbrExtraInteractiveManager, CbrProcessExit, MessageParser, MessageSerializer, MessageQueingPolicy, SendMessageResult, TMessageIn, TMessageOut
from . import CbrExtraInteractiveHandler as CbrExtraInteractiveHandlerModule

class CbrExtraInteractive(Generic[TMessageIn, TMessageOut]):

    def __init__(
        self,
        manager: CbrExtraInteractiveManager,
        parser: MessageParser[TMessageIn],
        serializer: MessageSerializer[TMessageOut],
        queing_policy: MessageQueingPolicy
    ):
        self.__manager: CbrExtraInteractiveManager = manager
        self.__handler: CbrExtraInteractiveHandler[TMessageIn, TMessageOut] = CbrExtraInteractiveHandler(manager, parser, serializer, queing_policy)
        self.__subscriptions = SubscriptionComposite([])

    def stop(self) -> None:
        self.__manager.stop()

    def __on_subscribe(self, s: SubscriptionBase):
        self.__subscriptions.append(s)

    def observe_values(self) -> Dsl[TMessageIn]:
        return self.__handler.observe_values().with_subscribe_callback(self.__on_subscribe)
    
    def observe_status(self) -> Dsl[CbrProcessExit]:
        return self.__manager.observe_status().with_subscribe_callback(self.__on_subscribe)
    
    def send_message(self, message: TMessageOut) -> SendMessageResult:
        return self.__handler.send_message(message)
    
    def dispose_subscriptions(self) -> None:
        self.__subscriptions.dispose()
    
def run_interactive(
    args : List[str],
    parser: MessageParser[TMessageIn],
    serializer: MessageSerializer[TMessageOut],
    queing_policy: MessageQueingPolicy,
    parent: Optional[QObject] = None
) -> CbrExtraInteractive[TMessageIn, TMessageOut]:

    return CbrExtraInteractive(
        CbrExtraInteractiveHandlerModule.run_interactive(args, parent),
        parser,
        serializer,
        queing_policy
    )