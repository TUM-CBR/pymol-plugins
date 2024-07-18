from PyQt5.QtCore import pyqtSignal, QObject
import pytest
import time
from typing import Callable, TypeVar

TResult = TypeVar('TResult')

class TestSupportObject(QObject):
    run_action = pyqtSignal(object)

    def __init__(self):
        super().__init__()
        self.run_action.connect(self.__on_run_action)

    def __on_run_action(self, action):
        action()

test_support_object = None

def __init_test_support_object():
    global test_support_object
    if test_support_object is None:
        test_support_object = TestSupportObject()

@pytest.fixture
def pymol_fixture():
    import pymol as pymol_module
    pymol_module.finish_launching()
    time.sleep(1)
    run_in_ui_thread(__init_test_support_object)
    from pmg_qt.pymol_qt_gui import window
    yield pymol_module
    test_support_object.run_action.emit(lambda: window.close())
    # run_in_ui_thread(lambda: window.close())

def run_in_ui_thread(fn: Callable[[], TResult]) -> TResult:
    from pymol import cmd
    print("abut")
    result = None
    def wrapper():
        nonlocal result
        print("running in UI")
        result = fn()

    cmd._call_in_gui_thread(wrapper)
    return result # pyright: ignore
