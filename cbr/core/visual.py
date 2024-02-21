from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import QComboBox, QPushButton
from pymol import cmd as pymol
from typing import Optional, Set

from ..core.pymol import objects, structure
from ..core.pymol.structure import StructureSelection

class StructureSelector(QObject):

    structure_changed = pyqtSignal(object)

    def __init__(
        self,
        structures_combo : QComboBox,
        refresh_button : QPushButton,
        copy_button : Optional[QPushButton] = None
    ):
        super(StructureSelector, self).__init__()
        self.__structures_combo = structures_combo
        self.__refresh_button = refresh_button
        self.__copy_button = copy_button
        
        refresh_button.clicked.connect(self.__on_refresh_clicked)

        if copy_button:
            copy_button.clicked.connect(self.__on_copy_clicked)

        self.__structures_combo.currentIndexChanged.connect(self.__on_item_selected)

    @pyqtSlot()
    def __on_item_selected(self):
        self.structure_changed.emit(self.currentSelection)

    @pyqtSlot()
    def __on_copy_clicked(self):
        selection = self.currentSelection

        if not selection:
            return

        sequence = structure.get_selection_sequece(selection.selection)
        clipboard = QApplication.clipboard()

        assert clipboard is not None, "This function must be used with a GUI"
        clipboard.setText(sequence)

    @pyqtSlot()
    def __on_refresh_clicked(self):
        self.refresh_structures()

    def refresh_structures(self):
        self.__structures_combo.clear()

        for (name, chain) in objects.iter_chains():

            segs: Set[str] = set()
            pymol.iterate(
                "model %s and chain %s" % (name, chain),
                'segs.add(segi)',
                space={'segs': segs}
            )

            for seg in segs:
                self.__structures_combo.addItem(
                    "%s/%s/%s" % (name, chain, seg),
                    StructureSelection(name, chain, seg)
                )

    @property
    def currentSelection(self) -> Optional[StructureSelection]:
        return self.__structures_combo.currentData()

def as_structure_selector(
    structures_combo : QComboBox,
    refresh_button : QPushButton,
    copy_button : Optional[QPushButton] = None
) -> StructureSelector:
    return StructureSelector(structures_combo, refresh_button, copy_button)