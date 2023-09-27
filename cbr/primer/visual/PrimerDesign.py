from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QFileDialog, QWidget
import re
from typing import Optional

from ...core.Context import Context
from ...core.Qt.QtCore import run_in_thread
from ...core.Qt.QtWidgets import progress_manager, show_error, show_exception
from .PrimerViewer import PrimerViewer
from .Ui_PrimerDesign import Ui_PrimerDesign

dna_re = re.compile(r"(?P<left>(A|C|T|G)+)\[(?P<design>(A|C|T|G)+)\](?P<right>(A|C|T|G)+)", re.IGNORECASE)

KNOWN_ORGANISMS = ['E_COLI', 'P_PASTORIS']

PrimerOrganism = str

class PrimerDesign(QWidget):

    def __init__(
        self,
        context : Context
    ):
        super(PrimerDesign, self).__init__()
        self.__ui = Ui_PrimerDesign()
        self.__ui.setupUi(self)
        self.__context = context

        self.__ui.designPrimersButton.clicked.connect(self.__on_design_primers)
        self.__progress = progress_manager(
            self.__ui.primerDesignProgress,
            self.__ui.designPrimersButton
        )
        self.__progress.on_result.connect(self.__on_design_primers_result)
        self.__progress.on_exception.connect(self.__on_design_primers_error)
        self.__ui.organismCombo.addItems(KNOWN_ORGANISMS)
        self.__ui.openPrimersButton.clicked.connect(self.__on_open_primers)

    def __get_organism(self) -> Optional[PrimerOrganism]:
        organism = self.__ui.organismCombo.currentText()

        if organism in KNOWN_ORGANISMS:
            return organism

        show_error(self, "Select a known organism, found: '%s'" % organism)

    @pyqtSlot()
    def __on_open_primers(self):
        result_file,_ = QFileDialog.getOpenFileName(
            self,
            "Select Primers Database File",
            "",
            f"Primers Database Files (*.sqlite)"
        )

        self.__on_design_primers_result(result_file)

    @pyqtSlot(object)
    def __on_design_primers_result(self, db: str):
        self.__context.run_widget(
            lambda _: PrimerViewer(db)
        ).show()

    @pyqtSlot(Exception)
    def __on_design_primers_error(self, exception : Exception):
        show_exception(self, exception)

    @pyqtSlot()
    def __on_design_primers(self):
        sequence = self.__ui.sequenceInputText.toPlainText()
        organism = self.__get_organism()
        if organism is None:
            return

        result = self.__design_primers(sequence, organism)
        self.__progress.watch_progress(result)
        
    @run_in_thread
    def __design_primers(
        self,
        raw_sequence : str,
        organism : PrimerOrganism
    ) -> str:

        sections = dna_re.match(raw_sequence)

        if sections is None:
            raise ValueError(invalid_sequence_error)

        left = sections.group('left')
        design = sections.group('design')
        right = sections.group('right')
        sequence = left + design + right

        raise Exception("Not Implemented!")

invalid_sequence_error =\
"""The sequence must only consists of C,T,G,A and must contain a region surrounded
by brackets ([...CTGA...]) which will be the region for which primers will be
designed.
"""