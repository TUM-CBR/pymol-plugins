from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt, pyqtSlot
from typing import Any, Iterable, List, Optional

from PyQt5.QtWidgets import QTableView

from ..QtWidgets import show_info

def all_the_keys(records : Iterable[dict]):
    return set(
        key
        for record in records
        for key in record.keys()
    )

class JsonRecordsModel(QAbstractTableModel):

    def __init__(self, data: List[dict]):
        super().__init__()
        self.__data = data
        self.__keys = list(all_the_keys(data))

    def records(self):
        return self.__data

    def rowCount(self, parent = None) -> int:
        return len(self.__data)

    def columnCount(self, parent = None) -> int:
        return len(self.__keys)

    def __get_headers(self) -> List[str]:
        return list(self.__keys)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            headers = self.__get_headers()
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def display_role_data(self, index : QModelIndex, truncate : Optional[int] = None):
        record = self.__data[index.row()]
        key = self.__keys[index.column()]
        value = record.get(key)

        if value is None:
            return None

        if isinstance(value, str):
            return value if truncate is None else value[:truncate]

        if isinstance(value, (int, bool, float)):
            return value

        return "{..}"

    def data(self, index: QModelIndex, role: int = Qt.DisplayRole) -> Any:

        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
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

    def set_records(self, records: List[dict]):
        self.__model = JsonRecordsModel(records)
        self.setModel(self.__model)

    def append_records(self, records: List[dict]):
        records = self.__model.records() + records
        self.__model = JsonRecordsModel(records)
        self.setModel(self.__model)