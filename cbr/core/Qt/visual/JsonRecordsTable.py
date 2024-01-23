from PyQt5.QtCore import QModelIndex, Qt, pyqtSlot
from typing import Any, Dict, Iterable, List, Optional

from PyQt5.QtWidgets import QTableView

from ..QtWidgets import show_info
from .value_handlers import EditRecordModelBase

Record = Dict[Any, Any]

def all_the_keys(records : Iterable[Record]):
    return set(
        key
        for record in records
        for key in record.keys()
    )

class JsonRecordsModel(EditRecordModelBase):

    def __init__(
        self,
        data: List[Record],
        read_only : bool = True
    ):
        super().__init__()
        self.__read_only = read_only
        self.set_data(data, notify=False)

    def set_data(self, data: List[Record], notify: bool = True):
        self.__data = data
        self.__keys = list(all_the_keys(data))

        if notify:
            self.modelReset.emit()

    def records(self):
        return self.__data

    def rowCount(self, parent : Any = None) -> int:
        return len(self.__data)

    def columnCount(self, parent : Any = None) -> int:
        return len(self.__keys)

    def __get_headers(self) -> List[str]:
        return list(self.__keys)

    def headerData(
        self,
        section: int,
        orientation : Qt.Orientation,
        role : int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            headers = self.__get_headers()
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)
    
    def __get_value__(self, index: QModelIndex) -> Optional[Any]:
        record = self.__data[index.row()]
        key = self.__keys[index.column()]
        return record.get(key)
    
    def __set_value__(self, index: QModelIndex, value: Any) -> bool:
        record = self.__data[index.row()]
        key = self.__keys[index.column()]

        record[key] = value
        return True

    def display_role_data(self, index : QModelIndex, truncate : Optional[int] = None):
        value = self.__get_value__(index)
        if value is None:
            return None

        if isinstance(value, str):
            return value if truncate is None else value[:truncate]

        if isinstance(value, (int, bool, float)):
            return value

        return "{..}"
    
    def flags(self, index: QModelIndex) -> Qt.ItemFlags:

        flags = super().flags(index)
        if self.__read_only:
            return flags

        if self.__can_edit__(index):
            return flags | Qt.ItemFlag.ItemIsEditable

        return flags

    def data(
        self,
        index: QModelIndex,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:

        if not index.isValid():
            return None

        if role == Qt.ItemDataRole.DisplayRole:
            return self.display_role_data(index, truncate=100)

class JsonRecordsTable(QTableView):

    def __init__(self):
        super().__init__()
        self.__model = JsonRecordsModel([])
        self.setModel(self.__model)
        self.doubleClicked.connect(self.__on_cell_double_clicked)

    @pyqtSlot(QModelIndex)
    def __on_cell_double_clicked(self, index : QModelIndex):
        data = self.__model.display_role_data(index)
        show_info(self, "Value", str(data))

    def set_records(self, records: List[Record]):
        self.__model = JsonRecordsModel(records)
        self.setModel(self.__model)

    def append_records(self, records: List[Record]):
        records = self.__model.records() + records
        self.__model = JsonRecordsModel(records)
        self.setModel(self.__model)