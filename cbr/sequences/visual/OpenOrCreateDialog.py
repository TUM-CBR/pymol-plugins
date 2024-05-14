from glob import glob
from os import path
from PyQt5.QtCore import Qt, pyqtSlot, QAbstractTableModel, QModelIndex
from PyQt5.QtWidgets import QDialog, QFileDialog, QWidget
import re
from typing import Any, List, NamedTuple, Optional, Sequence

from ...core.Qt.QtWidgets import show_error
from .Ui_CreateDatabaseDialog import Ui_CreateDatabaseDialog

FILENAME_RE = re.compile(r"(\w|\d)+")

class OpenOrCreateResult(NamedTuple):
    db_file: str
    search_folders: Sequence[str]

class DbListModel(QAbstractTableModel):

    K_COLUMN_NAME = "Database"
    COLUMNS = [K_COLUMN_NAME]

    def __init__(self, files: Sequence[str], parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.__files = files

    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__files)
    
    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.COLUMNS)
    
    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if role != Qt.ItemDataRole.DisplayRole:
            return None
        
        if orientation == Qt.Orientation.Horizontal:
            return self.COLUMNS[section]
        
    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if role == Qt.ItemDataRole.DisplayRole:
            return path.basename(self.__files[index.row()])
        
    def get_file(self, index: QModelIndex) -> str:
        return self.__files[index.row()]

class OpenOrCreateDialog(QDialog):

    def __init__(
        self,
        db_directory: str,
        parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(parent)
        self.__ui = Ui_CreateDatabaseDialog()
        self.__ui.setupUi(self)

        self.__search_folders: List[str] = []
        self.__db_directory: str = db_directory
        self.__result: Optional[OpenOrCreateResult] = None

        self.__ui.createButton.clicked.connect(self.__on_create_clicked)
        self.__ui.browseButton.clicked.connect(self.__select_search_folder)
        self.__ui.cancelButton.clicked.connect(self.__on_cancel)

        existing_files = glob(path.join(db_directory, "*.sqlite"))
        self.__db_list_model = DbListModel(existing_files)
        self.__ui.databasesTable.setModel(self.__db_list_model)

        self.__ui.databasesTable.doubleClicked.connect(self.__on_item_double_clicked)

    @pyqtSlot(QModelIndex)
    def __on_item_double_clicked(self, item: QModelIndex):
        self.__result = OpenOrCreateResult(
            self.__db_list_model.get_file(item),
            []
        )
        self.accept()

    def exec(self) -> int:
        self.__result = None
        self.__search_folders = []
        return super().exec()
    
    @pyqtSlot()
    def __select_search_folder(self):
        folder = QFileDialog.getExistingDirectory(self, caption="Select search folder")

        if path.exists(folder):
            self.__search_folders = [folder]
            self.__ui.selectedFileLabel.setText(path.basename(folder))

    @pyqtSlot()
    def __on_cancel(self):
        self.__result = None
        self.__search_folders = []
        self.reject()

    @pyqtSlot()
    def __on_create_clicked(self):

        name = self.__ui.nameLineEdit.text()
        if not FILENAME_RE.fullmatch(name):
            show_error(self, "New Database Error", f"Only letters and numbers are allowed in the database name.")
            return

        db_file = path.join(self.__db_directory, f"{name}.sqlite")
        if path.exists(db_file):
            show_error(self, "New Database Error", f"The database {name} already exits!")
            return
        
        if len(self.__search_folders) < 1:
            show_error(self, "Search Folder Error", f"You must select one search folder to create the database.")
            return

        self.__result = OpenOrCreateResult(db_file, self.__search_folders)

        self.accept()

    def value(self):
        return self.__result