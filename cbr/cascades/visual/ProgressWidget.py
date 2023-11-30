from typing import Optional
from PyQt5.QtCore import QProcess, pyqtSlot
from PyQt5.QtWidgets import QWidget

from ...core.Qt.visual.JsonRecordsTable import JsonRecordsTable
from ...extra.main import CbrExtraProcess
from .Ui_ProgressWidget import Ui_ProgressWidget

class ProgressWidget(QWidget):

    def __init__(self):
        super().__init__()
        self.__ui = Ui_ProgressWidget()
        self.__ui.setupUi(self)
        self.__process : Optional[CbrExtraProcess] = None

        self.__progress_details = JsonRecordsTable()
        self.layout().replaceWidget(
            self.__ui.detailsTable,
            self.__progress_details
        )
        self.__progress_details.setVisible(False)
        
        self.__ui.detailsButton.clicked.connect(self.__toggle_details)
        self.setVisible(False)

    @pyqtSlot()
    def __toggle_details(self):

        if self.__progress_details.isVisible():
            self.__ui.detailsButton.setText("Show Details")
            self.__progress_details.setVisible(False)
        else:
            self.__ui.detailsButton.setText("Hide Details")
            self.__progress_details.setVisible(True)

    def monitor_progress(self, process : CbrExtraProcess):

        self.setVisible(True)
        self.__progress_details.set_records([])
        self.__process = process

        self.__process.message_signal.connect(self.__on_message)
        self.__process.finished.connect(self.__on_finish)

    @pyqtSlot(int, QProcess.ExitStatus)
    def __on_finish(self, exitCode, exitStatus):
        self.setVisible(False)

    @pyqtSlot(object)
    def __on_message(self, message: dict):
        self.__progress_details.append_records([message])