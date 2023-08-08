import itertools
from typing import Iterable, List, NamedTuple, Tuple
from PyQt5.QtCore import QPoint, Qt, pyqtSlot
from PyQt5.QtWidgets import QHeaderView, QTableWidgetItem, QWidget

from ...core.Qt.QtWidgets import open_copy_context_menu
from .Ui_ChimerasResults import Ui_ChimerasResults

BASE_PAIR_COMPLEMENTS = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

def get_complement(base_pair : str) -> str:
    try:
        return BASE_PAIR_COMPLEMENTS[base_pair.upper()]
    except KeyError:
        raise ValueError("The character %s is not a base pair." % base_pair)

class Chimera(NamedTuple):
    name : str
    fragments : List[str]

class ChimeraGeneratorPosition(NamedTuple):
    sequences : List[str]

class ChimerasGeneratorArgs(NamedTuple):
    positions : List[ChimeraGeneratorPosition]
    binding_codons : List[str]
    overhang_length : int

    def __clean_binding_codons(self, sequence : str) -> str:
        for overhang in self.binding_codons:
            sequence = sequence.replace(overhang, "")

        return sequence

    def __clean_sequence(self, sequence : str) -> str:
        return self.__clean_binding_codons(sequence.upper())

    def __combine_sequences(self, seq_a : List[str], seq_b : str) -> List[str]:
        last = seq_a[-1]

        if(self.overhang_length > 0 and last[-self.overhang_length:] != seq_b[0:self.overhang_length]):
            raise ValueError("The sequence '%s' and '%s' have mismatched overhangs" % (last, seq_b))

        return seq_a + [seq_b[self.overhang_length:]]

    def __make_chimera(self, sequences : List[Tuple[int, str]]):
        chimera = "".join(str(i) for i,_ in sequences)

        result = [sequences[0][1]]
        for (_, sequence) in sequences[1:]:
            result = self.__combine_sequences(result, sequence)

        return Chimera(
            name = chimera,
            fragments = result
        )

    def __generate_chimeras(self):
        catalogue = [
            [(i + 1,self.__clean_sequence(sequence)) for i,sequence in enumerate(position.sequences)]
            for position in self.positions
        ]
        return [
            self.__make_chimera(list(sequences))
            for sequences in itertools.product(*catalogue)
        ]

    @staticmethod
    def generate_chimeras(
        positions : List[ChimeraGeneratorPosition],
        binding_codons : Iterable[str],
        overhang_length : int
    ) -> List[Chimera]:
        return ChimerasGeneratorArgs(
            positions = positions,
            binding_codons = list(binding_codons),
            overhang_length = overhang_length
        ).__generate_chimeras()

K_COMBINED = "Combined"
K_FRAGMENTS = "Fragments"

class ChimerasResults(QWidget):

    def __init__(self, chimeras : List[Chimera]):

        super(ChimerasResults, self).__init__()

        self.__chimeras = chimeras
        self.__ui = Ui_ChimerasResults()
        self.__ui.setupUi(self)
        self.__ui.resultsTable.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        self.__ui.resultsTable.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.__ui.resultsTable.customContextMenuRequested.connect(self.__on_results_context_menu)

        self.__ui.displayCombo.addItems([K_COMBINED, K_FRAGMENTS])
        self.__ui.displayCombo.currentTextChanged.connect(self.__update_chimeras)

        self.__render_chimeras()

    @pyqtSlot(QPoint)
    def __on_results_context_menu(self, pos):
        open_copy_context_menu(self.__ui.resultsTable, pos)

    @pyqtSlot(str)
    def __update_chimeras(self, mode : str):

        if mode == K_COMBINED:
            self.__render_chimeras()
        else:
            self.__render_fragments()

    def __render_fragments(self):
        resultsTable = self.__ui.resultsTable
        resultsTable.clear()
        column_count = 1 + len(self.__chimeras[0].fragments)
        resultsTable.setColumnCount(column_count)
        headers = ["Chimera"] + list(str(i) for i in range(1, column_count))
        resultsTable.setHorizontalHeaderLabels(headers)

        for i,chimera in enumerate(self.__chimeras):
            resultsTable.setItem(i, 0, QTableWidgetItem(chimera.name))
            for j,fragment in enumerate(chimera.fragments):
                resultsTable.setItem(i, j + 1, QTableWidgetItem(fragment))
    
    def __render_chimeras(self):

        resultsTable = self.__ui.resultsTable
        resultsTable.clear()
        resultsTable.setColumnCount(2)
        headers = ["Chimera", "Sequence"]
        resultsTable.setHorizontalHeaderLabels(headers)
        resultsTable.setRowCount(len(self.__chimeras))

        for i,chimera in enumerate(self.__chimeras):
            resultsTable.setItem(i, 0, QTableWidgetItem(chimera.name))
            resultsTable.setItem(i, 1, QTableWidgetItem("".join(chimera.fragments)))
