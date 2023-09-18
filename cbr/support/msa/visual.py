from PyQt5.QtCore import QObject, pyqtSlot
from PyQt5.QtWidgets import QComboBox, QPushButton
from typing import Optional, List, Tuple

from ...clustal.Clustal import Clustal
from ...clustal import msa
from ...core.pymol import structure
from ...core.pymol.structure import StructureSelection
from ...core.visual import as_structure_selector, StructureSelector

from .MsaSelector import Msa

class AlignmentIndexSelector(QObject):

    def __init__(
        self,
        structure_selector : StructureSelector,
        sequence_combo : QComboBox,
        msa : Msa,
        clustal : Clustal
    ):
        super(AlignmentIndexSelector, self).__init__()

        self.__sequence_combo = sequence_combo
        self.__sequence_combo.currentIndexChanged.connect(self.__on_sequence_changed)
        
        self.__structures_selector = structure_selector
        self.__structures_selector.structure_changed.connect(self.__on_structure_changed)

        self.__clustal = clustal
        self.__msa = msa
        self.__msa_relative_positions : Optional[List[int]] = None
        
        self.__update_msa_combo(msa)

    @pyqtSlot()
    def __on_structure_changed(self, structure : StructureSelection):
        self.__on_any_change()

    @pyqtSlot(int)
    def __on_sequence_changed(self, index : int):
        self.__on_any_change()

    @property
    def __msa_item(self) -> Optional[Tuple[str, str]]:
        return self.__sequence_combo.currentData()

    @property
    def relative_positions(self) -> Optional[List[int]]:
        return self.__msa_relative_positions

    def __on_any_change(self):

        structure_selection = self.__structures_selector.currentSelection

        if structure_selection is None:
            return

        structure_sequence = structure.get_pdb_sequence(structure_selection)

        msa_item = self.__msa_item

        if msa_item is None:
            return

        msa_sequence_name, msa_sequence = msa_item

        join_msa = self.__clustal.run_msa_items(
            [ ("structure", structure_sequence)
            , (msa_sequence_name, msa.clean_msa_blanks(msa_sequence))
            ]
        )

        self.__msa_relative_positions = list(msa.get_relative_positions(self.__msa, join_msa))

    def __update_msa_combo(self, msa : Msa):

        self.__sequence_combo.clear()

        for name, sequence in msa.items():
            self.__sequence_combo.addItem(
                name,
                (name, sequence)
            )

def as_alignment_link_selector(
    structures_combo : QComboBox,
    refresh_button : QPushButton,
    sequence_combo : QComboBox,
    msa : Msa,
    clustal : Clustal
):
    return AlignmentIndexSelector(
        as_structure_selector(structures_combo, refresh_button),
        sequence_combo,
        msa,
        clustal
    )