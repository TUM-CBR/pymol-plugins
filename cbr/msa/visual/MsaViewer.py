from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget
from typing import Dict, Optional

from ...clustal import Clustal
from ...core.Context import Context
from ...core import visual
from ...core.pymol.structure import StructureSelection
from ...support.msa import msa_selector, Msa

from .ColorByMotif import ColorByMotif
from .ColorByResidue import ColorByResidue
from .Ui_MsaViewer import Ui_MsaViewer
from .MsaContext import MsaContext

class MsaViewer(QWidget):

    NoStructureSelectedException = ValueError("No structure has been selected for the alignment")

    class MsaContextImpl(MsaContext):

        def __init__(
            self,
            msa_viewer : 'MsaViewer'
        ) -> None:
            super().__init__()
            self.__msa_viewer = msa_viewer

        @property
        def selected_msa_sequence_name(self) -> str:
            return self.__msa_viewer.__ui__.sequenceCombo.currentText()

        @property
        def selected_structure(self) -> Optional[StructureSelection]:
            return self.__msa_viewer.__structure_selector__.currentSelection

        @property
        def sequences(self) -> Msa:
            return self.__msa_viewer.__sequences__

    def __init__(self, context : Context, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.__ui = Ui_MsaViewer()
        self.__ui.setupUi(self)
        self.__structure_selector = visual.as_structure_selector(
            self.__ui.structuresCombo,
            self.__ui.structuresRefreshButton
        )
        self.__set_sequences({})
        self.__clustal = Clustal.get_clustal_from_context(context)
        self.__msa_selector = msa_selector(
            self,
            self.__ui.selectMsaButton,
            self.__ui.selectedFileLabel,
            self.__ui.createMsaButton
        )
        self.__msa_selector.msa_file_selected.connect(self.__on_select_msa)

        msa_context = self.MsaContextImpl(self)
        self.__ui.msaColoringWidgets.addTab(
            ColorByResidue(self.__clustal, msa_context),
            "Color by Residue"
        )

        self.__ui.msaColoringWidgets.addTab(
            ColorByMotif(msa_context),
            "Color by Motif"
        )

    @property
    def __structure_selector__(self) -> visual.StructureSelector:
        return self.__structure_selector

    @property
    def __ui__(self) -> Ui_MsaViewer:
        return self.__ui

    @property
    def __sequences__(self) -> Msa:
        return self.__sequences

    @pyqtSlot(object)
    def __on_select_msa(self, _):
        self.__set_sequences(self.__msa_selector.msa or {})

    def __set_sequences(self, sequences: Msa):
        self.__sequences = sequences
        self.__ui.sequenceCombo.clear()
        self.__ui.sequenceCombo.addItems(self.__sequences.keys())

