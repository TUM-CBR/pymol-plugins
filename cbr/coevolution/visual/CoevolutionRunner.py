from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QFileDialog, QWidget
from typing import Optional

from cbr.coevolution.visual.CoevolutionViewer import CoevolultionViewer

from ...acpsicov import main as acpsicov
from ...core.Context import Context
from ...core.Qt.QtCore import run_in_thread
from ...core.Qt.QtWidgets import progress_manager, show_error, show_exception
from ...support import msa

from .Ui_CoevolutionRunner import Ui_CoevolutionRunner

ACPSICOV_EXTENSIONS = "*.acpsicov"

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
        self.__progress_manager.on_exception.connect(self.__on_exception)
        self.__context = context
        self.__ui.openButton.clicked.connect(self.__on_open_result)

    @pyqtSlot()
    def __on_open_result(self):

        msa = self.__msa_selector.msa
        if msa is None:
            show_error(self, "You must select an MSA file first.")
            return

        result_file,_ = QFileDialog.getOpenFileName(
            self,
            "Select ACPSICOV File",
            "",
            f"ACPSICOV files ({ACPSICOV_EXTENSIONS})"
        )

        try:
            with open(result_file, "r") as acpsicov_file:
                result = acpsicov.acpsicov_load(acpsicov_file)
                self.__context.run_widget(
                    lambda ctx: CoevolultionViewer(
                        result,
                        msa,
                        ctx
                    )
                ).show()

        except Exception as e:
            show_exception(self, e)

    @pyqtSlot()
    def __on_exception(self, exn : Exception):
        show_exception(self, exn)

    @pyqtSlot()
    def __on_run_coevolution(self):

        self.__progress_manager.watch_progress(self.__run_coevolution(self.__msa_selector.msa))
        pass

    @run_in_thread
    def __run_coevolution(self, msa : Optional[msa.Msa]):

        if msa is None:
            raise Exception("Please select an MSA file!")

        acpsicov.acpsicov(msa)