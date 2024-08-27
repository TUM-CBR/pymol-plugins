from pymol import cmd
from PyQt5.QtCore import pyqtSlot, QPoint, Qt
from PyQt5.QtWidgets import QTableWidgetItem, QWidget
from typing import Callable, Dict, List, NamedTuple, Optional, Set, Tuple

from ...clustal import msa
from ...clustal.Clustal import Clustal
from ...core.pymol.structure import StructureSelection
from ...core.pymol.visual.PymolResidueResultsTable import PymolResidueSelector
from ...core.Qt.QtWidgets import open_copy_context_menu, with_error_handler
from ...core.sequence import residue_to_3
from ...support.msa import Msa

from .MsaContext import MsaContext
from .support import msa_to_structure_position_map, sequence_to_structure_position_map
from .Ui_ColorByResidue import Ui_ColorByResidue

COLOR_MAX = 999
COLOR_MIN = 600

def perc_str(n : float):
    return str(round(100*n, ndigits=2))

K_GAP = "-"

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
        if K_GAP in score:
            del score[K_GAP]
        self.__structure_residue = structure_residue
        self.__score = score
        self.__pdb_position = pdb_position
        self.__msa_position = msa_position

        self.__score_perc = self.__compute_score_percent()

    @property
    def msa_consensus(self) -> str:
        scores = self.score_percent
        rank = list(scores.keys())
        rank.sort(key=lambda x: scores[x], reverse=True)
        best = rank[0]
        score = scores[best]
        best_score = f"({perc_str(score)}%)"

        if score < self.blanks_percent:
            best = residue_to_3(rank[1]) + "*"
        else:
            best = residue_to_3(best)

        return f"{best} {best_score}"

    @property
    def blanks(self) -> int:
        return self.__blanks

    @property
    def blanks_percent(self) -> float:
        return self.blanks / self.sum_with_blanks

    @property
    def residue(self) -> str:
        resi = self.__structure_residue
        score = self.__score_perc.get(resi) or 0.0
        resi = residue_to_3(resi)
        return f"{resi} ({perc_str(score)}%)"

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
        return self.__score_perc

    def __compute_score_percent(self) -> Dict[str, float]:
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

class ColorByResidue(QWidget):

    scoring_methods : Dict[str, Callable[[MsaConservationResult], float]] = {
        'maximum conservation': lambda x: x.max_conservation_score,
        'maximum conservation (gaps weightd)': lambda x: x.max_conservation_weighted,
        'residue count': lambda x: x.residue_count_score,
        'relative to structure': lambda x: x.relative_to_structure_conservation_score,
        'structure relative to consensus': lambda x: x.structure_relative_to_consensus_score
    }

    def __init__(
        self,
        clustal : Clustal,
        msa_context : MsaContext
    ) -> None:
        super().__init__()

        self.__msa_context = msa_context
        self.__clustal = clustal

        self.__ui = Ui_ColorByResidue()
        self.__ui.setupUi(self)

        self.__ui.scoringCombo.addItems(ColorByResidue.scoring_methods)

        self.__residue_selector = PymolResidueSelector(self.__ui.resultsTable)

        # Context menu for copy/paste
        self.__ui.resultsTable.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.__ui.resultsTable.customContextMenuRequested.connect(self.__on_results_context_menu)
        self.__ui.colorConservedButton.clicked.connect(self.__on_colorConservedButton_clicked)

    @property
    def __selected_structure(self) -> Optional[StructureSelection]:
        return self.__msa_context.selected_structure

    @property
    def __sequences(self) -> Msa:
        return self.__msa_context.sequences

    @property
    def __selected_structure_name(self) -> str:
        return self.__msa_context.selected_msa_sequence_name

    def __get_structure_sequence(self) -> str:
        return self.__msa_context.get_structure_sequence()

    @pyqtSlot(QPoint)
    def __on_results_context_menu(self, pos):
        open_copy_context_menu(self.__ui.resultsTable, pos)

    def __get_structure_positions(self)  -> List[int]:
        return msa_to_structure_position_map(
            self.__selected_structure_name,
            self.__sequences,
            self.__get_structure_sequence()
        )

    def __get_structure_conservation(self) -> MsaConservationResults:
        selection = self.__selected_structure
        if selection is None:
            raise Exception("You must select a structure to use this feature!")

        sequences = list(self.__sequences.values())
        positions = self.__get_structure_positions()
        structure_sequence = self.__get_structure_sequence()
        result = {}
        results_keys = set()
        offset = sequence_to_structure_position_map(selection)

        for (i, ix) in enumerate(positions):

            if ix in result:
                continue

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

            result[ix] = MsaConservationResult(offset[ix], i+1, structure_residue, conserved)

        return MsaConservationResults(results_keys, result)
    
    @property
    def __scoring_function(self) -> Callable[[MsaConservationResult], float]:
        return self.scoring_methods.get(self.__ui.scoringCombo.currentText()) or (lambda x: 0.0)

    def __get_columns(self, residues : Set[str]):

        def get_name(res : str) -> str:
            res = res.upper()
            if res == 'X':
                return "Skipped 'X'"
            elif res == "-":
                return "Gap"

            return residue_to_3(res)

        return list(sorted(((key, get_name(key)) for key in residues), key=lambda kv: kv[1]))

    def __set_results_table(self, results : MsaConservationResults) -> None:
        table = self.__ui.resultsTable
        table.reset()
        residue_selector = self.__residue_selector
        residue_selector.reset()
        scoring = self.__scoring_function
        structure_selection = self.__selected_structure

        if structure_selection is None:
            raise Exception("A structure must be selected to perform this action!")

        self.__residue_selector.set_selection(structure_selection.selection)

        base_columns = ['structure position', 'MSA position', 'structure residue', 'consensus residue', 'score', 'Gaps (%)']

        columns = self.__get_columns(results.results_keys)
        column_headers = [col[1] for col in columns]

        column_names = base_columns + \
            ["%s (%%)" % r for r in column_headers]
        table.setColumnCount(len(column_names))
        table.setRowCount(len(results.results))
        table.setHorizontalHeaderLabels(column_names)

        for (i,score) in enumerate(results.results.values()):
            perc_scores = score.score_percent
            consensus = score.msa_consensus

            residue_selector.set_item_riesidues(i, [score.pdb_position])
            table.setItem(i, 0, QTableWidgetItem(str(score.pdb_position)))
            table.setItem(i, 1, QTableWidgetItem(str(score.msa_position)))
            table.setItem(i, 2, QTableWidgetItem(score.residue))
            table.setItem(i, 3, QTableWidgetItem(consensus))
            table.setItem(i, 4, QTableWidgetItem(str(scoring(score))))

            table.setItem(
                i,
                5,
                QTableWidgetItem("%s%%" % perc_str(score.blanks_percent))
            )

            for (j, (key, _)) in enumerate(columns):
                k_score = perc_str(perc_scores.get(key) or 0.0)
                table.setItem(i, j + len(base_columns), QTableWidgetItem("%s%%" % k_score))


    @pyqtSlot(int)
    def on_gapTresholdSlider_sliderMoved(self, value : int):
        self.__ui.gapTresholdLabel.setText("%i%%" % value)

    @with_error_handler()
    def __color_by_conservation(self):
        selection = self.__selected_structure

        if selection is None:
            raise Exception("A structure must be selected to perform this action!")

        conservation = self.__get_structure_conservation()
        scoring = self.__scoring_function
        scores = [(result, scoring(result)) for (i, result) in conservation.results.items()]
        factor = (COLOR_MAX - COLOR_MIN) / max(score for (r, score) in scores)
        structure_query = selection.selection
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
    def __on_colorConservedButton_clicked(self):
        self.__color_by_conservation()
