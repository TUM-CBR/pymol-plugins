import os
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QDoubleValidator, QIntValidator
from PyQt5.QtWidgets import QFileDialog, QWidget
import random
import re
from tempfile import TemporaryDirectory
from typing import Optional

from ...core.Context import Context
from ...core.Qt.QtCore import run_in_thread
from ...core.Qt.QtWidgets import progress_manager, show_error, show_exception
from ..operations import design_primers
from .PrimerViewer import PrimerViewer
from .Ui_PrimerDesign import Ui_PrimerDesign

dna_re = re.compile(r"(?P<left>(A|C|T|G)+)\[(?P<design>(A|C|T|G)+)\](?P<right>(A|C|T|G)+)", re.IGNORECASE)

KNOWN_ORGANISMS = ['E_COLI', 'P_PASTORIS']

PrimerOrganism = str

DEFAULT_PRIMER_ARGS = {

}

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

        # Primer size validators
        validator = QIntValidator(1,50)
        self.__ui.minSizeInput.setValidator(validator)
        self.__ui.maxSizeInput.setValidator(validator)

        # DNA concentration validator
        validator = QDoubleValidator(0.01, 10000, 2)
        self.__ui.concInput.setValidator(validator)

        self.__progress = progress_manager(
            self.__ui.primerDesignProgress,
            self.__ui.designPrimersButton
        )
        self.__progress.on_result.connect(self.__on_design_primers_result)
        self.__progress.on_exception.connect(self.__on_design_primers_error)
        self.__ui.organismCombo.addItems(KNOWN_ORGANISMS)
        self.__ui.openPrimersButton.clicked.connect(self.__on_open_primers)
        self.__workdir = TemporaryDirectory()

    def __del__(self):
        self.__workdir.cleanup()

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

        self.__open_primers_db(result_file)

    @pyqtSlot(object)
    def __on_design_primers_result(self, db: str):
        self.__open_primers_db(db, delete_when_closed=True)

    def __open_primers_db(self, db: str, delete_when_closed: bool = False):
        self.__context.run_widget(
            lambda _: PrimerViewer(db, delete_when_closed=delete_when_closed)
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

        min_size = int(self.__ui.minSizeInput.text())
        max_size = int(self.__ui.maxSizeInput.text())

        if min_size >= max_size:
            show_error(self, "Input Error", "Minimum primer size must be smaller than maximum primer size.")

        result = self.__design_primers(sequence, organism, min_size, max_size)
        self.__progress.watch_progress(result)

    def __new_file_name(self):
        dir = self.__workdir.name
        return os.path.join(
            dir,
            f"primers_{len(os.listdir(dir))}_{random.randint(0,2**16)}.sqlite"
        )
        
    @run_in_thread
    def __design_primers(
        self,
        raw_sequence : str,
        organism : PrimerOrganism,
        min_primer_size,
        max_primer_size
    ) -> str:

        sections = dna_re.match(raw_sequence)

        if sections is None:
            raise ValueError(invalid_sequence_error)

        left = sections.group('left')
        design = sections.group('design')
        right = sections.group('right')
        sequence = left + design + right
        result_file = self.__new_file_name()
        design_primers(
            len(left),
            int(len(design)/3),
            min_primer_size,
            max_primer_size,
            organism,
            sequence,
            result_file
        )
        return result_file

invalid_sequence_error =\
"""The sequence must only consists of C,T,G,A and must contain a region surrounded
by brackets ([...CTGA...]) which will be the region for which primers will be
designed.
"""