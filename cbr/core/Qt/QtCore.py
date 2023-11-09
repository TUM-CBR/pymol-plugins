from concurrent.futures import Future
from PyQt5.QtCore import QAbstractTableModel, QObject, pyqtSlot, QModelIndex, Qt, QTimer, QThread, QObject, pyqtSignal, pyqtSlot
import os
from typing import Any, Callable, Dict, List, TypeVar

class Throttle(QObject):

    def __init__(
        self,
        timeout : int,
        action,
        *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.__timer = QTimer()
        self.__timer.setSingleShot(True)
        self.__timer.timeout.connect(self.__on_timeout)
        self.__timeout = timeout
        self.__action = action
        self.__args = []
        self.__kwargs = {}

    @pyqtSlot()
    def __on_timeout(self):
        self.__action(*self.__args, **self.__kwargs)

    def trigger(self, *args, **kwargs):
        self.__args = args
        self.__kwargs = kwargs
        self.__timer.start(self.__timeout)

TResult = TypeVar("TResult")

def debug_qthread():
    if "QTHREAD_DEBUGGING" in os.environ:
        import debugpy # pyright: ignore
        debugpy.debug_this_thread()

def run_in_thread(func: Callable[..., TResult]) -> Callable[..., 'Future[TResult]']:
    def wrapper(me, *args: List[Any], **kwargs: dict) -> Future:
        future = Future()

        class Worker(QThread):

            finished = pyqtSignal()

            def run(self) -> None:
                debug_qthread()
                try:
                    result = func(me, *args, **kwargs)
                    future.set_result(result)
                except Exception as e:
                    future.set_exception(e)
                self.finished.emit()

        thread = Worker(me)
        thread.finished.connect(thread.deleteLater)

        #thread.started.connect(worker.run)
        thread.start()

        return future

    return wrapper

class DictionaryModel(QAbstractTableModel):

    def __init__(
        self,
        values : Dict[str, str]
    ):
        super().__init__()
        self.__values = list(values.items())

    def rowCount(self, parent = None) -> int:
        return len(self.__values)

    def columnCount(self, parent = None) -> int:
        return 2

    def data(self, index: QModelIndex, role=Qt.DisplayRole):

        if not index.isValid():
            return None

        row_ix = index.row()
        col_ix = index.column()
        item = self.__values[row_ix]

        if role == Qt.DisplayRole:
            return item[col_ix]


class BlockSignalsContextManager(object):

    def __init__(self, obj : QObject):
        self.__obj = obj

    def __enter__(self, *args, **kwargs):
        self.__obj.blockSignals(True)
        return self

    def __exit__(self, *args, **kwargs):
        self.__obj.blockSignals(False)

block_signals = BlockSignalsContextManager