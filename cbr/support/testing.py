from PyQt5.QtCore import pyqtSignal, pyqtSlot, QObject
from PyQt5.QtWidgets import QMainWindow
import time
from typing import Any, Callable, Optional, TypeVar

TResult = TypeVar('TResult')

class RunActionRequest:
    action: Callable[[], Any]
    result: Optional[Any] = None
    error: Optional[Exception] = None
    done: bool = False

    def __init__(self, action: Callable[[], Any]):
        self.action = action

    def wait(self) -> Any:
        while not self.done:
            time.sleep(0.1)
        if self.error is not None:
            raise self.error
        return self.result

class UiRunner(QObject):

    singleton: Optional['UiRunner'] = None

    __run_action = pyqtSignal(RunActionRequest)

    __close_request = pyqtSignal()

    def __init__(self):
        super().__init__()
        from pmg_qt.pymol_qt_gui import window
        self.__window: Optional[QMainWindow] = window
        self.__close_request.connect(self.__on_close_request)
        self.__run_action.connect(self.__on_run_action)

    @pyqtSlot(RunActionRequest)
    def __on_run_action(self, request: RunActionRequest):
        import threading
        try:
            print("UI thread: ", threading.current_thread().ident)
            request.result = request.action()
        except Exception as e:
            request.error = e
        finally:
            request.done = True

    @classmethod
    def init_singleton(cls):
        print("init singleton")

        assert cls.singleton is None, "Only one singleton can be created"

        cls.singleton = cls()

    def run_in_ui_async(self, action: Callable[[], TResult]) -> RunActionRequest:

        if self.__window is None:
            raise Exception("The Ui runner cannot be used after the application has been closed.")

        request = RunActionRequest(action)
        self.__run_action.emit(request)
        return request

    def __on_close_request(self):
        window = self.__window

        if window is None:
            return
        self.__window = None
        window.close()

    def close(self):
        self.__close_request.emit()

class PymolTestSupportObject(QObject):
 
    def __init__(self):
        super().__init__()

        while UiRunner.singleton is None:
            time.sleep(0.1)
            print("waiting for the singleton")
        self.__runner = UiRunner.singleton

    def run_in_ui_async(self, action: Callable[[], TResult]) -> RunActionRequest:
        return self.__runner.run_in_ui_async(action)

    def run_in_ui(self, action: Callable[[], TResult]) -> TResult:
        request = self.run_in_ui_async(action)
        return request.wait()

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self.__runner.close()
        import time
        time.sleep(1)


