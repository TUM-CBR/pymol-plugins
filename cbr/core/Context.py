from tempfile import TemporaryDirectory
from PyQt5.QtWidgets import QMainWindow, QWidget
from PyQt5.QtGui import QCloseEvent
from typing import Any, Callable, Dict, List, TypeVar

from .executable import *

T = TypeVar('T')

class RunQWidgetContext(QMainWindow):
    
    def __init__(
        self,
        widget : QWidget,
        on_close : 'Callable[[RunQWidgetContext], None]'
    ):
        super().__init__()

        self.__on_close = on_close
        self.setCentralWidget(widget)

    def closeEvent(self, a0: QCloseEvent) -> None:
        super().closeEvent(a0)
        self.__on_close(self)

class Context(object):

    def __init__(self):
        self.__store : Dict[str, Any] = {}
        self.__widgets : List[RunQWidgetContext] = []
        self.__directories : List[TemporaryDirectory[Any]] = []

    def create_temporary_directory(self) -> str:
        directory = TemporaryDirectory()
        self.__directories.append(directory)
        return directory.name
    
    def __cleanup(self):
        for directory in self.__directories:
            directory.cleanup()

        self.__directories = []

    def __del__(self):
        self.__cleanup()

    def get_executable(
        self,
        widget: QWidget,
        executable: KnownExecutables
    ) -> Executable:
        
        import shutil
        if shutil.which("protein_mpnn_run.sh"):
            return Executable(
                executable_type=ExecutableType.Binary,
                location="protein_mpnn_run.sh"
            )
        
        raise Exception("Cannot locate the exectuable!")
                

    def create_or_load(self, name : str, load : Callable[[], T]) -> T:

        if name in self.__store:
            return self.__store[name]

        result = self.__store[name] = load()

        return result

    def __on_close(self, window: RunQWidgetContext):

        try:
            self.__widgets.remove(window)
        except ValueError:
            # This should not happen, but not critical if it does
            pass

    def run_widget(self, run_widget : Callable[['Context'], QWidget]) -> RunQWidgetContext:

        widget = run_widget(self)
        window = RunQWidgetContext(widget, self.__on_close)

        self.__widgets.append(window)

        return window