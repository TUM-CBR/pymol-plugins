from os import path
from PyQt5.QtCore import pyqtSlot, QRegularExpression
from PyQt5.QtWidgets import QDialog, QFileDialog, QWidget
from PyQt5.QtGui import QRegExpValidator
from typing import List, NamedTuple, Optional, Sequence

from ...core.Qt.QtWidgets import show_error
from .Ui_CreateDatabaseDialog import Ui_CreateDatabaseDialog

FILENAME_RE = r"(\w\d)+"

class OpenOrCreateResult(NamedTuple):
    db_file: str
    search_folders: Sequence[str]

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
        self.__ui.nameLineEdit.setValidator(QRegExpValidator(QRegularExpression(FILENAME_RE)))
        self.__result: Optional[OpenOrCreateResult] = None

        self.__ui.createButton.clicked.connect(self.__on_create_clicked)
        self.__ui.browseButton.clicked.connect(self.__select_search_folder)
        self.__ui.cancelButton.clicked.connect(self.__on_cancel)

    def exec(self) -> int:
        self.__result = None
        self.__search_folders = []
        return super().exec()
    
    @pyqtSlot()
    def __select_search_folder(self):
        folder = QFileDialog.getExistingDirectory(self, caption="Select search folder")

        if path.exists(folder):
            self.__search_folders = [folder]

    @pyqtSlot()
    def __on_cancel(self):
        self.__result = None
        self.__search_folders = []
        self.reject()

    @pyqtSlot()
    def __on_create_clicked(self):

        name = self.__ui.nameLineEdit.text()
        db_file = path.join(self.__db_directory, f"{name}.sqlite")

        if path.exists(db_file):
            show_error(self, "New Database Error", f"The database {name} already exits!")
            return
        
        if len(self.__search_folders):
            show_error(self, "Search Folder Error", f"You must select one search folder to create the database.")

        self.__result = OpenOrCreateResult(db_file, self.__search_folders)

        self.accept()

    def value(self):
        return self.__result