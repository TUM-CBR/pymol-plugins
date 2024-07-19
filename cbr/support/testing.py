from PyQt5.QtCore import pyqtSignal, pyqtSlot, QObject
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

    run_action = pyqtSignal(RunActionRequest)

    def __init__(self):
        super().__init__()
        self.run_action.connect(self.__on_run_action)
        self.__valid = True

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
        cls.singleton = cls()

    def invalidate(self):
        self.__valid = False

    def run_in_ui_async(self, action: Callable[[], TResult]) -> RunActionRequest:

        if not self.__valid:
            raise Exception("The Ui runner cannot be used after the application has been closed.")

        request = RunActionRequest(action)
        self.run_action.emit(request)
        return request

class PymolTestSupportObject(QObject):
 
    def __init__(self):
        super().__init__()

        self.__window = None

        while self.__window is None:
            from pmg_qt.pymol_qt_gui import window
            time.sleep(0.1)
            self.__window = window

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

    def __cleanup_ui(self):
        print("closing window")
        self.__window.close()
        print('goodbye')

    def __exit__(self, *args, **kwargs):
        self.run_in_ui_async(self.__cleanup_ui)
        self.__runner.invalidate()
        import time
        time.sleep(1)


