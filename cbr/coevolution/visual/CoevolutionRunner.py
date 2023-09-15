from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget
from typing import Optional

from ...acpsicov import main as acpsicov
from ...core.Context import Context
from ...core.Qt.QtCore import run_in_thread
from ...core.Qt.QtWidgets import progress_manager, show_error, show_exception
from ...support import msa

from .Ui_CoevolutionRunner import Ui_CoevolutionRunner

class CoevolutionRunner(QWidget):

    def __init__(self, context : Context) -> None:
        super().__init__()
        self.__ui = Ui_CoevolutionRunner()
        self.__ui.setupUi(self)
        self.__msa_selector = msa.msa_selector(
            self,
            self.__ui.selectFileButton,
            self.__ui.selectedFileLabel
        )

        self.__ui.runCoevolutionButton.clicked.connect(self.__on_run_coevolution)
        self.__progress_manager = progress_manager(
            self.__ui.coevolutionProgress,
            self.__ui.runCoevolutionButton
        )

    @pyqtSlot()
    def __on_run_coevolution(self):

        self.__progress_manager.watch_progress(self.__run_coevolution(self.__msa_selector.msa))
        pass

    @run_in_thread
    def __run_coevolution(self, msa : Optional[msa.Msa]):

        if msa is None:
            raise Exception("Please select an MSA file!")

        acpsicov.acpsicov(msa)