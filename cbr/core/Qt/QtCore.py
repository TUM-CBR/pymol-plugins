from concurrent.futures import Future
from PyQt5.QtCore import QAbstractTableModel, QObject, pyqtSlot, QModelIndex, Qt, QTimer, QThread, QObject, pyqtSignal, pyqtSlot
import os
from typing import Any, Callable, Dict, Generic, List, TypeVar

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

    def __init__(self, *obj : QObject):
        self.__obj = obj

    def __enter__(self, *args, **kwargs):

        for obj in self.__obj:
            obj.blockSignals(True)
        return self

    def __exit__(self, *args, **kwargs):

        for obj in self.__obj:
            obj.blockSignals(False)

block_signals = BlockSignalsContextManager

TKey = TypeVar("TKey")
TModel = TypeVar("TModel")

class AbstractRecordTableModel(QAbstractTableModel, Generic[TKey, TModel]):

    def __init__(
        self,
        values: List[TModel]
        ) -> None:
        super().__init__()

        kv = [(self.get_uid(v), v) for v in values]
        self.__keys = [k for k, _v in kv]
        self.__values = dict(kv)

    def get_record(self, ix: int) -> TModel:
        return self.__values[self.__keys[ix]]

    def get_record_count(self) -> int:
        return len(self.__keys)

    def add_records(self, *records: TModel):
        for record in records:
            uid = self.get_uid(record)

            if uid not in self.__values:
                self.__keys.append(uid)
            
            self.__values[uid] = record

    def clear(self) -> None:
        self.beginResetModel()
        self.__keys = []
        self.__values = {}
        self.endResetModel()

    def get_uid(self, value : TModel) -> TKey:
        raise NotImplementedError("RecordModel requires an implementation of get_uid")
