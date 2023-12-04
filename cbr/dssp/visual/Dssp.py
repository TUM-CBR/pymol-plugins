from datetime import datetime
from os import path
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget
import tempfile

from ...core.Qt.QtWidgets import show_error, with_error_handler
from ...core.visual import as_structure_selector
from .Ui_dssp import Ui_Dssp

class Dssp(QWidget):

    def __init__(self) -> None:
        super().__init__()

        self.__ui = Ui_Dssp()
        self.__ui.setupUi(self)

        self.__structure_selector = as_structure_selector(
            self.__ui.structuresCombo,
            self.__ui.refreshButton
        )

        self.__ui.runButton.clicked.connect(self.__on_run)

        self.__working_directory = tempfile.TemporaryDirectory()

    def __del__(self):
        self.__working_directory.cleanup()

    @pyqtSlot(name="__on_run")
    @with_error_handler()
    def __on_run(self):

        name = datetime.now().strftime('%Y_%m_%d_%H%M%S')
        pdb_file = path.join(
            self.__working_directory.name,
            f"{name}.pdb"
        )

        selection = self.__structure_selector.currentSelection

        if selection is None:
            show_error(
                self,
                "Error",
                "No structure has been selected!"
            )
            return

        
