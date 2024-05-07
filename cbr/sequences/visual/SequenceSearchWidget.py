from glob import glob
from os import path
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QDialog, QWidget
from typing import Optional, Sequence

from ...core.Context import Context
from ...core import namedtuple
from ..data import SequencesConfig
from .OpenOrCreateDialog import OpenOrCreateDialog
from .Ui_SequenceSearchWidget import Ui_SequenceSearchWidget

class SequenceSearchWidget(QWidget):

    def __init__(
        self,
        context: Context,
        parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(parent)
        self.__context = context
        self.__config_directory = context.get_config_directory("sequences")
        self.__config_file = path.join(self.__config_directory, "settings.json")
        self.__ui = Ui_SequenceSearchWidget()
        self.__ui.setupUi(self)
        self.__set_busy(False)

        self.__database: Optional[str] = None
        self.__open_dialog = OpenOrCreateDialog(self.__databases_directory())
        self.__open_default_database()

    def __set_busy(self, busy: bool):
        self.__ui.scanProgress.setVisible(busy)
        self.__ui.rescanButton.setEnabled(not busy)
        self.__ui.searchButton.setEnabled(not busy)
        self.__ui.changeDatabase.setEnabled(not busy)

    def __databases_directory(self):
        return self.__config_directory
    
    def __list_databases(self) -> Sequence[str]:
        pattern = path.join(self.__databases_directory(), '*.sqlite')
        return list(glob(pattern))
    
    def __read_config(self) -> SequencesConfig:

        try:
            with open(self.__config_file, 'r') as config:
                return namedtuple.load(SequencesConfig, config)
        except FileNotFoundError:
            return SequencesConfig()
        except ValueError:
            return SequencesConfig()

    def __scan(self, search_folders: Optional[Sequence[str]] = None):
        return

    def __open_or_create_database(self) -> bool:

        if self.__open_dialog.exec() == QDialog.DialogCode.Accepted:
            value = self.__open_dialog.value()
            assert value is not None, "Bug in OpenOrCreateDialog, returned None upon success."
            self.__database = value.db_file
            self.__scan(value.search_folders)
            return True
        else:
            return False

    def __open_default_database(self):

        config = self.__read_config()

        if config.selected_database is not None:
            self.__database = config.selected_database
        elif not self.__open_or_create_database():
            self.__context.close_app()