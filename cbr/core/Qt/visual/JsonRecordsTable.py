from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt
from typing import Any, Iterable, List

from PyQt5.QtWidgets import QTableView

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

    def __display_role_data(self, index : QModelIndex):
        record = self.__data[index.row()]
        key = self.__keys[index.column()]
        value = record.get(key)

        if value is None:
            return None

        if isinstance(value, (str, int, bool, float)):
            return value

        return "{..}"

    def data(self, index: QModelIndex, role: int = Qt.DisplayRole) -> Any:

        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            return self.__display_role_data(index)

class JsonRecordsTable(QTableView):

    def __init__(self):
        super().__init__()
        self.__model = JsonRecordsModel([])
        self.setModel(self.__model)

    def set_records(self, records: List[dict]):
        self.__model = JsonRecordsModel(records)
        self.setModel(self.__model)

    def append_records(self, records: List[dict]):
        records = self.__model.records() + records
        self.__model = JsonRecordsModel(records)
        self.setModel(self.__model)