from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget
import re
from typing import List, NamedTuple, Optional

from ...core.Context import Context
from ...core.Qt.QtCore import run_in_thread
from ...core.Qt.QtWidgets import progress_manager, show_error, show_exception
from ..Operations import CODONS_MAP, DesignPrimersResult, Operations, PrimerOrganism
from .PrimerViewer import PrimerViewer
from .Ui_PrimerDesign import Ui_PrimerDesign

dna_re = re.compile(r"(?P<left>(A|C|T|G)+)\[(?P<design>(A|C|T|G)+)\](?P<right>(A|C|T|G)+)", re.IGNORECASE)

class PrimersDesign(NamedTuple):
    primers : List[DesignPrimersResult]
    sequence : str

class PrimerDesign(QWidget):

    def __init__(
        self,
        context : Context,
        operations : Optional[Operations] = None
    ):
        super(PrimerDesign, self).__init__()
        self.__ui = Ui_PrimerDesign()
        self.__ui.setupUi(self)
        self.__operations = operations or Operations()
        self.__context = context

        self.__ui.designPrimersButton.clicked.connect(self.__on_design_primers)
        self.__progress = progress_manager(
            self.__ui.primerDesignProgress,
            self.__ui.designPrimersButton
        )
        self.__progress.on_result.connect(self.__on_design_primers_result)
        self.__progress.on_exception.connect(self.__on_design_primers_error)
        self.__ui.organismCombo.addItems(CODONS_MAP.keys())

    def __get_organism(self) -> Optional[PrimerOrganism]:
        organism = self.__ui.organismCombo.currentText()

        if organism in CODONS_MAP:
            return organism

        show_error(self, "Select a known organism, found: '%s'" % organism)

    @pyqtSlot(object)
    def __on_design_primers_result(self, result : PrimersDesign):
        print(result)
        self.__context.run_widget(
            lambda _: PrimerViewer(result.primers, result.sequence)
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
    ) -> PrimersDesign:

        sections = dna_re.match(raw_sequence)

        if sections is None:
            raise ValueError(invalid_sequence_error)

        left = sections.group('left')
        design = sections.group('design')
        right = sections.group('right')
        sequence = left + design + right
        results = self.__operations.design_primers(
            sequence,
            len(left),
            len(design),
            organism
        )
        return PrimersDesign(
            primers=list(results),
            sequence=sequence
        )

invalid_sequence_error =\
"""The sequence must only consists of C,T,G,A and must contain a region surrounded
by brackets ([...CTGA...]) which will be the region for which primers will be
designed.
"""