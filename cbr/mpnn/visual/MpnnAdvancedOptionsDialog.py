from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QDialog

from ...core.Qt.visual.NamedTupleEditor import namedtuple_eidtor
from ..data import MpnnArgs
from .Ui_MpnnAdvancedOptions import Ui_MpnnAdvancedOptions

class MpnnAdvancedOptionsDialog(QDialog):

    def __init__(self):
        super().__init__()

        self.__ui = Ui_MpnnAdvancedOptions()
        self.__ui.setupUi(self)

        self.__options_editor = namedtuple_eidtor(
            self.__ui.mpnnArgsTable,
            MpnnArgs()
        )

        self.__ui.resetButton.clicked.connect(self.__on_reset_clicked)
        self.__ui.acceptButton.clicked.connect(self.__on_accept_clicked)

    @pyqtSlot()
    def __on_reset_clicked(self):
        self.__options_editor[0] = MpnnArgs()

    @pyqtSlot()
    def __on_accept_clicked(self):
        self.accept()

    def value(self) -> MpnnArgs:
        result = self.__options_editor[0]
        assert result is not None, "Advanced options shold never be None"
        return result