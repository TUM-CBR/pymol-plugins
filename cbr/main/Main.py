from . import resources

from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget

from ..core.Context import Context
from ..msa import Applications as Msa
from ..schema import Applications as Schema
from .Ui_Main import Ui_Main

class Main(QWidget):

    def __init__(self, context : Context, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        self.__context = context
        self.__ui = Ui_Main()
        self.__ui.setupUi(self)

    @pyqtSlot()
    def on_schemaRasppButton_clicked(self):
        print("starting raspp")
        Schema.schema_raspp(self.__context)

    @pyqtSlot()
    def on_msaViewerButton_clicked(self):
        print("starting msa viewer")
        Msa.msa_viewer(self.__context)
