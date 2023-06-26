from abc import ABC, abstractmethod, abstractproperty
from typing import Callable, Generic, TypeVar


class IContext(ABC):

    @abstractproperty
    def task_manager(self) -> 'ITaskManager':
        ...

TaskValue = TypeVar("TaskValue")
TaskSpec = Callable[[], TaskValue]

class ITaskInstance(ABC, Generic[TaskValue]):

    @abstractproperty
    def error(self) -> Exception:
        ...

    @abstractproperty
    def result(self) -> 'TaskValue | None':
        ...

    @abstractproperty
    def is_working(self) -> bool:
        ...

    @abstractmethod
    def on_completed(self, action : Callable[[], None]) -> None:
        pass

class ITaskManager(ABC):

    @abstractmethod
    def run_task(self, name : str, task_spec : TaskSpec[TaskValue]) -> ITaskInstance[TaskValue]:
        pass