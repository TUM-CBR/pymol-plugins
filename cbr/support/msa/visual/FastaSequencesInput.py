from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QDialog, QWidget
from typing import Dict

from ....clustal import Clustal
from .Ui_FastaSequencesInput import Ui_FastaSequencesInput

class FastaSequencesInput(QDialog):

    def __init__(self, parent : QWidget):
        super(FastaSequencesInput, self).__init__(parent=parent)
        self.__clustal = Clustal.get_clustal()
        self.__msa_result : 'Dict[str,str] | None' = None
        self.__ui = Ui_FastaSequencesInput()
        self.__ui.setupUi(self)

    @property
    def msa_result(self) -> 'Dict[str,str] | None':
        return self.__msa_result

    @pyqtSlot()
    def on_sequenceButtonBox_accepted(self):

        try:
            self.__msa_result = self.__clustal.run_msa_items(self.__ui.sequenceTextEdit.toPlainText())
            self.__ui.errorLabel.setText(str())
            self.done(QDialog.Accepted)
        except ValueError as error:
            self.__ui.errorLabel.setText(str(error))

    @pyqtSlot()
    def on_sequenceButtonBox_rejected(self):
        self.done(QDialog.Rejected)
    

    