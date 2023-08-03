import pymol.cmd as cmd
from PyQt5.QtCore import pyqtSlot, QPoint, Qt
from PyQt5.QtWidgets import QDialog, QFileDialog, QTableWidgetItem, QWidget
from typing import Callable, Dict, List, NamedTuple, Set, Tuple

from ...clustal import msa
from ...clustal import Clustal
from ...core.Context import Context
from ...core import visual
from ...core.pymol import selection
from ...core.pymol import structure
from ...core.pymol.visual.PymolResidueResultsTable import PymolResidueSelector
from ...core.Qt.QtWidgets import open_copy_context_menu

from .FastaSequencesInput import FastaSequencesInput
from .Ui_MsaViewer import Ui_MsaViewer

COLOR_MAX = 999
COLOR_MIN = 600

def perc_str(n : float):
    return str(round(100*n, ndigits=2))

class MsaConservationResult(object):

    def __init__(
            self,
            pdb_position : int,
            msa_position : int,
            structure_residue : str,
            score : Dict[str, int]
        ):
        score = score.copy()
        self.__blanks = score.get("-") or 0
        if "-" in score:
            del score["-"]
        self.__structure_residue = structure_residue
        self.__score = score
        self.__pdb_position = pdb_position
        self.__msa_position = msa_position

    @property
    def blanks(self) -> int:
        return self.__blanks

    @property
    def blanks_percent(self) -> float:
        return self.blanks / self.sum_with_blanks

    @property
    def residue(self) -> str:
        return self.__structure_residue

    @property
    def pdb_position(self) -> int:
        return self.__pdb_position

    @property
    def msa_position(self) -> int:
        return self.__msa_position

    @property
    def score(self) -> Dict[str, int]:
        return self.__score
    
    @property
    def sum_with_blanks(self) -> float:
        return sum(self.__score.values()) + self.__blanks

    @property
    def score_percent(self) -> Dict[str, float]:
        total = self.sum_with_blanks
        return dict([
            (res, score/total)
            for (res, score) in self.__score.items()
        ])

    @property
    def residue_count_score(self) -> float:
        return len(self.__score)
    
    @property
    def max_conservation_score(self) -> int:
        return max(self.__score.values())
    
    @property
    def max_conservation_weighted(self) -> int:
        return self.max_conservation_score + self.__blanks
    
    @property
    def relative_to_structure_conservation_score(self) -> int:
        return self.__score[self.__structure_residue]

    @property
    def structure_relative_to_consensus_score(self):
        return self.relative_to_structure_conservation_score / self.max_conservation_score
    
class MsaConservationResults(NamedTuple):
    results_keys : Set[str]
    results : Dict[int, MsaConservationResult]

class MsaViewer(QWidget):

    scoring_methods : Dict[str, Callable[[MsaConservationResult], float]] = {
        'maximum conservation': lambda x: x.max_conservation_score,
        'maximum conservation (gaps weightd)': lambda x: x.max_conservation_weighted,
        'residue count': lambda x: x.residue_count_score,
        'relative to structure': lambda x: x.relative_to_structure_conservation_score,
        'structure relative to consensus': lambda x: x.structure_relative_to_consensus_score
    }

    def __init__(self, context : Context, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.__ui = Ui_MsaViewer()
        self.__ui.setupUi(self)
        self.__ui.scoringCombo.addItems(MsaViewer.scoring_methods)
        visual.as_structure_selector(self.__ui.structuresCombo, self.__ui.structuresRefreshButton)
        self.__set_sequences({})
        self.__clustal = Clustal.get_clustal_from_context(context)
        self.__msa_input_dialog = FastaSequencesInput(context, parent = self)
        self.__residue_selector = PymolResidueSelector(self.__ui.resultsTable)

        # Context menu for copy/paste
        self.__ui.resultsTable.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.__ui.resultsTable.customContextMenuRequested.connect(self.__on_results_context_menu)

    @pyqtSlot(QPoint)
    def __on_results_context_menu(self, pos):
        open_copy_context_menu(self.__ui.resultsTable, pos)

    def __set_sequences(self, sequences : Dict[str, str]):
        self.__sequences = dict(sequences)
        self.__ui.sequenceCombo.clear()
        self.__ui.sequenceCombo.addItems(self.__sequences.keys())

    def on_createMsaButton_clicked(self):

        if self.__msa_input_dialog.exec() == QDialog.Accepted:
            assert self.__msa_input_dialog.msa_result, "Bug in the code, result should be set"
            self.__set_sequences(self.__msa_input_dialog.msa_result)

    @pyqtSlot()
    def on_selectMsaButton_clicked(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly  # Optionally set file dialog options

        # Show the file dialog and get the selected file
        file_path, _ = QFileDialog.getOpenFileName(self, "Open File", "", "FASTA alignment Files (*.*)", options=options)

        if not file_path:
            return

        print("Selected file:", file_path)
        self.__set_sequences(msa.parse_alignments(file_path))

    @property
    def __selected_structure(self) -> Tuple[str, str]:
        return self.__ui.structuresCombo.currentData()

    def __get_structure_sequence(self):
        (structure_name, chain) = self.__selected_structure
        return structure.get_pdb_sequence(structure_name, chain)

    def __get_structure_positions(self)  -> 'List[int | None]':
        sequence_name = self.__ui.sequenceCombo.currentText()
        sequence = msa.clean_msa_blanks(self.__sequences[sequence_name])
        structure_sequence = self.__get_structure_sequence()
        result = self.__clustal.run_msa_items(
            [ ("structure", structure_sequence)
            , (sequence_name, sequence)
            ]
        )

        return list(msa.get_relative_positions(self.__sequences, result))

    def __get_structure_conservation(self) -> MsaConservationResults:
        sequences = list(self.__sequences.values())
        positions = self.__get_structure_positions()
        structure_sequence = self.__get_structure_sequence()
        result = {}
        results_keys = set()
        (selected_structure, chain) = self.__selected_structure
        offset = structure.get_structure_offset(selected_structure, chain)

        for (i, ix) in filter(lambda x: x[1], enumerate(positions)):

            assert ix, "Bug in code, only positions that occur in the structure should be there"
            structure_residue = structure_sequence[ix].lower()
            conserved = {structure_residue: 1}
            for sequence in sequences:
                aa = sequence[i].lower()

                if not(msa.is_blank(aa)):
                    results_keys.add(aa)

                if aa in conserved:
                    conserved[aa] += 1
                else:
                    conserved[aa] = 1

            assert ix, "Bug in filtering criteria!"
            result[ix] = MsaConservationResult(ix + offset, i, structure_residue, conserved)

        return MsaConservationResults(results_keys, result)
    
    @property
    def __scoring_function(self) -> Callable[[MsaConservationResult], float]:
        return MsaViewer.scoring_methods.get(self.__ui.scoringCombo.currentText()) or (lambda x: 0.0)

    def __set_results_table(self, results : MsaConservationResults) -> None:
        table = self.__ui.resultsTable
        table.reset()
        residue_selector = self.__residue_selector
        residue_selector.reset()
        scoring = self.__scoring_function
        structure_selection = selection.get_selection_for_model(*self.__selected_structure)
        self.__residue_selector.set_selection(structure_selection)

        table.setColumnCount(len(results.results_keys) + 4)
        table.setRowCount(len(results.results))
        table.setHorizontalHeaderLabels(
            ['structure position', 'MSA position', 'structure residue', 'score'] + \
            ["%s (%%)" % r for r in results.results_keys] + \
            ["gaps (%)"]
        )

        for (i,score) in enumerate(results.results.values()):
            residue_selector.set_item_riesidues(i, [score.pdb_position])
            table.setItem(i, 0, QTableWidgetItem(str(score.pdb_position)))
            table.setItem(i, 1, QTableWidgetItem(str(score.msa_position)))
            table.setItem(i, 2, QTableWidgetItem(score.residue))
            table.setItem(i, 3, QTableWidgetItem(str(scoring(score))))

            for (j, key) in enumerate(results.results_keys):
                k_score = perc_str(score.score_percent.get(key) or 0.0)
                table.setItem(i, j + 4, QTableWidgetItem("%s%%" % k_score))

            table.setItem(
                i,
                len(results.results_keys) + 3,
                QTableWidgetItem("%s%%" % perc_str(score.blanks_percent))
            )

    @pyqtSlot(int)
    def on_gapTresholdSlider_sliderMoved(self, value : int):
        self.__ui.gapTresholdLabel.setText("%i%%" % value)

    def __color_by_conservation(self):
        (selected_structure, chain) = self.__selected_structure
        conservation = self.__get_structure_conservation()
        scoring = self.__scoring_function
        scores = [(result, scoring(result)) for (i, result) in conservation.results.items()]
        factor = (COLOR_MAX - COLOR_MIN) / max(score for (r, score) in scores)
        structure_query = structure.get_structure_query(selected_structure, chain)
        gaps_treshold = self.__ui.gapTresholdSlider.value()

        cmd.color("999", structure_query)

        for (result, score) in scores:
            ix = result.pdb_position
            if result.blanks_percent > gaps_treshold / 100:
                color = "white"
            else:
                i_score = int(score * factor)
                color = str(999 - i_score)
            cmd.color(color, "%s and resi %i-%i" % (structure_query, ix, ix))

        self.__set_results_table(conservation)

    @pyqtSlot()
    def on_colorConservedButton_clicked(self):
        self.__color_by_conservation()
