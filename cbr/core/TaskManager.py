from typing import Any, Callable, Generic, TypeVar
from PyQt5.QtCore import QObject, QThread, pyqtSignal, pyqtSlot

from .Context import Context

TaskValue = TypeVar("TaskValue")
TaskSpec = Callable[[], Any]

class TaskInstance(QObject, Generic[TaskValue]):

    work_started = pyqtSignal()
    work_completed = pyqtSignal()

    def __init__(
        self,
        task_id : int,
        task_name : str,
        task : TaskSpec,
        thread : QThread,
        *args,
        **kwargs
    ) -> None:
        super().__init__(*args, **kwargs)
        self.__task = task
        self.__is_working = False
        self.__task_id = task_id
        self.__result = None
        self.__error = None

    @property
    def is_working(self):
        return self.__is_working

    @property
    def error(self):
        return self.__error

    @property
    def result(self) -> 'TaskValue | None':
        return self.__result

    @pyqtSlot()
    def do_work(self, task_id : int):

        if(task_id != self.__task_id):
            return

        self.__is_working = True
        self.work_started.emit()
        pass

        try:
            self.__task()
        except Exception as e:
            self.__error = e
            print(e)

        self.work_completed.emit()

class TaskManager(QObject):

    start_work = pyqtSignal(int, name="start_work")

    def __init__(self, context : Context, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.__worker_thread = QThread()
        self.__worker_thread.start()
        self.__current_task_id = 0

    def __del__(self):
        self.__worker_thread.quit()

    def __next_task_id(self):
        self.__current_task_id += 1
        return self.__current_task_id

    def run_task(self, name : str, task_spec : TaskSpec) -> TaskInstance:
        tid = self.__next_task_id()
        task = TaskInstance(tid, name, task_spec, self.__worker_thread)
        task.moveToThread(self.__worker_thread)
        self.start_work.connect(task.do_work)
        self.start_work.emit(tid)

        return task

