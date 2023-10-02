from io import StringIO
import os
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt, pyqtSlot
from PyQt5.QtGui import QDoubleValidator, QIntValidator
from PyQt5.QtWidgets import QFileDialog, QWidget
import random
import re
from tempfile import TemporaryDirectory
from typing import Iterable, Optional

from ...core.Context import Context
from ...core.Qt.QtCore import run_in_thread
from ...core.Qt.QtWidgets import progress_manager, show_error, show_exception
from ...core.Qt.visual.EditRecordsDialog import EditRecordsDialog

from ..data import DEFAULT_PRIMER3_ARGS, Primer3Args
from ..operations import design_primers
from .PrimerViewer import PrimerViewer
from .Ui_PrimerDesign import Ui_PrimerDesign

dna_re = re.compile(r"(?P<left>(A|C|T|G)+)\[(?P<design>(A|C|T|G)+)\](?P<right>(A|C|T|G)+)", re.IGNORECASE)

KNOWN_ORGANISMS = ['E_COLI', 'P_PASTORIS']
GENE_SEQUENCE_DATA_ROLE = Qt.UserRole

PrimerOrganism = str

START_CODONS = ["ATG"]

STOP_CODONS = ["TAA", "TAG", "TGA"]

class GeneDataModel(QAbstractTableModel):

    def __init__(self, genes: Iterable[str]):
        super().__init__()
        self.__genes = [gene for gene in genes if len(gene) > 0]
        self.__columns = ["Size", "Sequence"]

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            headers = self.__columns
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def rowCount(self, parent=None):
        return len(self.__genes)

    def columnCount(self, parent=None):
        return len(self.__columns)

    def data(self, index: QModelIndex, role=Qt.DisplayRole):

        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            gene = self.__genes[index.row()]
            cols = [
                str(len(gene)),
                gene
            ]
            return cols[index.column()]
        elif role == GENE_SEQUENCE_DATA_ROLE:
            return self.__genes[index.row()]
        else:
            return None

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
        self.__ui.minSizeInput.setText("5")
        self.__ui.minSizeInput.setValidator(validator)
        self.__ui.maxSizeInput.setText("30")
        self.__ui.maxSizeInput.setValidator(validator)

        # Automatically find genes
        self.__ui.sequenceInputText.textChanged.connect(self.__on_text_changed)

        # Primer3 Aargs
        self.__primer3_args : Primer3Args = DEFAULT_PRIMER3_ARGS
        self.__ui.advancedButton.clicked.connect(self.__on_advanced_clicked)
        #self.__advanced_dialog = EditRecordsDialog(self, self.__primer3_args, DEFAULT_PRIMER3_ARGS)
        # DNA concentration validator
        validator = QDoubleValidator(0.01, 10000, 2)
        self.__ui.concInput.setValidator(validator)
        self.__ui.concInput.setText(str(self.__primer3_args.dna_conc))

        self.__progress = progress_manager(
            self.__ui.primerDesignProgress,
            self.__ui.designPrimersButton
        )
        self.__progress.on_result.connect(self.__on_design_primers_result)
        self.__progress.on_exception.connect(self.__on_design_primers_error)
        self.__ui.organismCombo.addItems(KNOWN_ORGANISMS)
        self.__ui.openPrimersButton.clicked.connect(self.__on_open_primers)
        self.__workdir = TemporaryDirectory()

    def __find_genes(self, plasmid: str) -> Iterable[str]:

        plasmid = plasmid.upper()
        pos = 3
        step = 1

        def next_codon():
            nonlocal pos
            nonlocal step
            nonlocal plasmid
            codon = plasmid[pos-3:pos]
            pos += step
            return codon

        while(pos <= len(plasmid)):

            codon = next_codon()
            if codon not in START_CODONS:
                continue

            step = 3

            # Beforehand we were reading one by one, so
            # next_codon had advanced by one on the last call.
            # We need to advance two more positions to be aligned
            # with the first codon of the gene
            pos += 2
            codon = next_codon()
            with StringIO() as gene:

                while(codon not in STOP_CODONS and pos <= len(plasmid)):
                    gene.write(codon)
                    codon = next_codon()

                gene.seek(0)

                if codon in STOP_CODONS:
                    yield gene.read()
                step = 1


    @pyqtSlot()
    def __on_text_changed(self):

        genes = self.__find_genes(self.__ui.sequenceInputText.toPlainText())
        self.__ui.selectGeneTable.setModel(GeneDataModel(genes))
        self.__ui.selectGeneTable.resizeColumnsToContents()

    @pyqtSlot()
    def __on_advanced_clicked(self):
        dialog = EditRecordsDialog(self, self.__primer3_args, DEFAULT_PRIMER3_ARGS)
        #self.__advanced_dialog.exec_()
        dialog.exec_()
        self.__primer3_args = dialog.new_value()
        self.__ui.concInput.setText(str(self.__primer3_args.dna_conc))

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

    def __get_plasmid_sequence(self) -> str:
        sequence = self.__ui.sequenceInputText.toPlainText()

        if dna_re.match(sequence):
            return sequence

        selection = list(self.__ui.selectGeneTable.selectedIndexes())

        if len(selection) != 1:
            raise Exception("Only one gene can be selected! Make sure only one cell is selected.")

        selected = selection[0]

        if selected:
            gene = selected.data(GENE_SEQUENCE_DATA_ROLE)
            return sequence.upper().replace("[", "").replace("]", "").replace(gene, f"[{gene}]")

        return sequence

    @pyqtSlot()
    def __on_design_primers(self):

        try:
            sequence = self.__get_plasmid_sequence()
            organism = self.__get_organism()
            if organism is None:
                return

            min_size = int(self.__ui.minSizeInput.text())
            max_size = int(self.__ui.maxSizeInput.text())

            if min_size >= max_size:
                show_error(self, "Input Error", "Minimum primer size must be smaller than maximum primer size.")

            result = self.__design_primers(sequence, organism, min_size, max_size, self.__primer3_args)
            self.__progress.watch_progress(result)
        except Exception as e:
            show_exception(self, e)

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
        min_primer_size : int,
        max_primer_size : int,
        primer3_args : Primer3Args
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
            primer3_args,
            result_file
        )
        return result_file

invalid_sequence_error =\
"""The sequence must only consists of C,T,G,A and must contain a region surrounded
by brackets ([...CTGA...]) which will be the region for which primers will be
designed.
"""