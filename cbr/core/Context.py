from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QMainWindow, QWidget
from typing import Any, Callable, Dict, TypeVar

T = TypeVar('T')

class RunQWidgetContext(object):
    
    def __init__(
        self, window : QMainWindow,
        widget : QWidget):

        self.__window = window
        self.__widget = widget

    @property
    def window(self) -> QMainWindow:
        return self.__window

    def show(self) -> None:

        self.__window.show()

class Context(object):

    def __init__(self):
        self.__store : Dict[str, Any] = {}
        self.__widgets = []

    def create_or_load(self, name : str, load : Callable[[], T]) -> T:

        if name in self.__store:
            return self.__store[name]

        result = self.__store[name] = load()

        return result

    def run_widget(self, run_widget : Callable[['Context'], QWidget]) -> RunQWidgetContext:

        window = QMainWindow()
        widget = run_widget(self)
        window.setCentralWidget(widget)

        instance = RunQWidgetContext(window, widget)
        self.__widgets.append(instance)

        return instance