from os import path
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QFileDialog, QWidget

from ...core.Qt.QtWidgets import with_error_handler
from .Ui_FastaViewerApp import Ui_FastaViewerApp

class FastaViewerApp(QWidget):

    def __init__(self) -> None:
        super().__init__()

        self.__ui = Ui_FastaViewerApp()
        self.__ui.setupUi(self)

        self.__ui.browseButton.clicked.connect(self.__on_browse_button_clicked)

    @pyqtSlot(name="__on_browse_button_clicked")
    @with_error_handler()
    def __on_browse_button_clicked(self):

        filename,_ = QFileDialog.getOpenFileName(self, "Select Fasta", filter="Fasta Files (*.fa *.fasta)")
        self.__ui.selectedFileLabel.setText(path.basename(filename))
        self.__ui.fastaViewerWidget.set_fasta(filename)

