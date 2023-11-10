from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QDialog, QWidget
from typing import cast, Optional

from ....core.pymol.structure import StructureSelection
from ....core.visual import as_structure_selector

from .Ui_MsaStructureSelector import Ui_MsaStructureSelector

class MsaStructureSelector(QDialog):

    def __init__(self, parent : QWidget):
        super().__init__(parent)
        self.__ui = Ui_MsaStructureSelector()
        self.__ui.setupUi(self)

        self.__structure_selector = as_structure_selector(
            self.__ui.structureCombo,
            self.__ui.refreshButton
        )

        self.__ui.confirmButtons.accepted.connect(self.__on_accept)
        self.__ui.confirmButtons.rejected.connect(self.__on_reject)

    def current_selection(self) -> Optional[StructureSelection]:
        return self.__structure_selector.currentSelection

    @pyqtSlot()
    def __on_reject(self):
        self.__ui.structureCombo.setCurrentIndex(-1)
        self.reject()

    @pyqtSlot()
    def __on_accept(self):
        self.accept()

    def exec(self) -> QDialog.DialogCode:
        self.__structure_selector.refresh_structures()
        return cast(QDialog.DialogCode, super().exec())