from os import path
import os
from PyQt5.QtCore import QDir
from PyQt5.QtWidgets import QMainWindow, QWidget
from PyQt5.QtGui import QCloseEvent
from tempfile import TemporaryDirectory
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

    def closeEvent(self, a0: Optional[QCloseEvent]) -> None:
        super().closeEvent(a0)
        self.__on_close(self)

class Context(object):

    def __init__(self):
        self.__store : Dict[str, Any] = {}
        self.__widgets : List[RunQWidgetContext] = []
        self.__directories : List[TemporaryDirectory[Any]] = []
        self.__application_directory = path.join(
            QDir.homePath(),
            ".cbrtools"
        )

    def __ensure_directory_exists(self, directory: str) -> str:
        if path.exists(directory):
            return directory
        
        os.mkdir(directory)
        return directory

    def __get_application_directory(self):
        return self.__ensure_directory_exists(self.__application_directory)

    def __components_directory(self):
        return self.__ensure_directory_exists(
            path.join(self.__get_application_directory(), "components")
        )

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
        
        return find_executable(
            widget,
            executable,
            self.__components_directory()
        )
                

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