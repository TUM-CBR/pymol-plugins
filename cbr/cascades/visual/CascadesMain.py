from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QFileDialog, QWidget

from ...core.Context import Context

from .CascadesViewer import CascadesViewer
from .Ui_CascadesMain import Ui_CascadesMain

class CascadesMain(QWidget):

    def __init__(
        self,
        context : Context
    ):
        super().__init__()
        self.__ui = Ui_CascadesMain()
        self.__context = context
        self.__ui.setupUi(self)

        self.__ui.loadExistingButton.clicked.connect(self.__on_select_file)

    @pyqtSlot()
    def __on_select_file(self):
        result_file,_ = QFileDialog.getOpenFileName(
            self,
            "Select a cascade BLAST",
            "",
            "Cascade BLAST databases (*.sqlite)"
        )

        if result_file is not None:
            self.__open_db(result_file)

    def __open_db(self, db_file: str):
        self.__context.run_widget(
            lambda _: CascadesViewer(db_file)
        ).show()