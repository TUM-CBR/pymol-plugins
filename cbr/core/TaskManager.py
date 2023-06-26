from PyQt5.QtCore import QMutex, QObject, QThread, pyqtSignal, pyqtSlot
from typing import Callable, Generic, TypeVar

from .Context import Context

TaskValue = TypeVar("TaskValue")
TaskSpec = Callable[[], TaskValue]

class QMutexWrapper():

    def __init__(self):
        self.__mutex = QMutex()

    def __enter__(self, *args, **kwargs) -> 'QMutexWrapper':
        self.__mutex.lock()
        return self

    def __exit__(self, *args, **kwargs) -> None:
        self.__mutex.unlock()

class TaskInstance(QObject, Generic[TaskValue]):

    NOT_STARTED = 0
    STARTED = 1
    COMPLETED = 2

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
        self.__state = TaskInstance.NOT_STARTED
        self.__task_id = task_id
        self.__result = None
        self.__error = None
        self.__mutex = QMutexWrapper()

    @property
    def is_working(self):
        return self.__state == TaskInstance.STARTED

    @property
    def error(self):
        return self.__error

    @property
    def result(self) -> 'TaskValue | None':
        return self.__result

    def __set_task_running(self):
        with self.__mutex:
            self.__state = TaskInstance.STARTED

    def __set_task_completed(self):
        with self.__mutex:
            self.__state = TaskInstance.COMPLETED

    def on_started(self, action: Callable[[], None]) -> None:

        def handler():
            try:
                action()
            finally:
                self.work_started.disconnect(handler)

        with self.__mutex:
            if self.is_working:
                should_emit = True
            else:
                should_emit = False
                self.work_started.connect(handler)

        if should_emit:
            action()

    def on_completed(self, action : Callable[[], None]) -> None:
        """
        This function will call action exactly once. If the
        task has completed, it will be called imediately, otherwise
        it will be called once the task has completed.
        """

        def handler():
            try:
                action()
            finally:
                self.work_completed.disconnect(handler)

        with self.__mutex:
            if self.__state == TaskInstance.COMPLETED:
                should_emit = True
            else:
                should_emit = False
                self.work_completed.connect(handler)

        if should_emit:
            action()

    @pyqtSlot()
    def do_work(self):

        with self.__mutex:
            if(self.__state == TaskInstance.COMPLETED):
                return

        self.__set_task_running()
        self.work_started.emit()
        pass

        try:
            self.__task()
        except Exception as e:
            self.__error = e
            print(e)
        finally:
            self.__set_task_completed()
            self.work_completed.emit()

class TaskManager(QObject):

    start_work = pyqtSignal(int)

    def __init__(self, context : Context, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.__worker_thread = QThread()
        self.__worker_thread.start()
        self.__current_task_id = 0

    @staticmethod
    def from_context(context: Context) -> 'TaskManager':
        return context.create_or_load(
            'TaskManager',
            lambda: TaskManager(context)
            )

    def __del__(self):
        self.__worker_thread.quit()

    def __next_task_id(self):
        self.__current_task_id += 1
        return self.__current_task_id

    def run_task(self, name : str, task_spec : TaskSpec[TaskValue]) -> TaskInstance[TaskValue]:
        tid = self.__next_task_id()
        task = TaskInstance(tid, name, task_spec, self.__worker_thread)
        task.moveToThread(self.__worker_thread)
        self.start_work.connect(task.do_work)
        self.start_work.emit(tid)
        task.on_completed(lambda: self.start_work.disconnect(task.do_work))

        return task

