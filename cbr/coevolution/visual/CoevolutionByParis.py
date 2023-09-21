from PyQt5.QtCore import QModelIndex, pyqtSlot
from PyQt5.QtGui import QDoubleValidator
from typing import Dict, Iterable, List, NamedTuple, Optional, Tuple

from PyQt5.QtWidgets import QTableWidget, QTableWidgetItem

from ...acpsicov.main import AcpsicovEntry, AcpsicovResult
from ...core.pymol.structure import ResidueDistanceCalculator, StructureSelection
from ...core.pymol import structure
from ...core.pymol.visual.PymolResidueResultsTable import pymol_residues_color_selector, ResidueWithIntensity
from ...support.msa import enumerate_pairs, Msa

from .CoevolutionResultViewerBase import CoevolutionResultViewerBase
from .Ui_CoevolutionByPairs import Ui_CoevolutionByPairs

class VisualizationStructure(NamedTuple):
    structure : StructureSelection
    link : List[int]
    sequence_index : Dict[int, str]

    @property
    def sequence_offset(self):
        return min(self.sequence_index.keys())

    def msa_to_structure_position(self, msa_pos: int) -> int:
        return self.link[msa_pos - 1]

    def msa_to_structure_pymol_posiiton(self, msa_pos: int) -> int:
        return self.msa_to_structure_position(msa_pos) + self.sequence_offset

    def msa_to_structure_residue(self, msa_pos: int) -> str:
        return self.sequence_index[self.msa_to_structure_pymol_posiiton(msa_pos)]

class ResultItemModel(NamedTuple):
    pymol_position_1 : Optional[int]
    pymol_position_2 : Optional[int]
    entry : AcpsicovEntry
    distance : Optional[float]

    @property
    def pos_1_fallback(self) -> int:
        return self.pymol_position_1 if self.pymol_position_1 is not None else self.entry.position1

    @property
    def pos_2_fallback(self) -> int:
        return self.pymol_position_2 if self.pymol_position_2 is not None else self.entry.position2

    @property
    def distance_str(self) -> str:
        return str(self.distance) if self.distance is not None else NA_VALUE

RESULT_ROLE = 1
ResultRoleType = ResultItemModel
NA_VALUE = "N/A"

class CoevolutionByPairs(CoevolutionResultViewerBase):

    def __init__(
        self,
        result : AcpsicovResult,
        msa : Msa
    ):

        super(CoevolutionByPairs, self).__init__()
        self.__ui = Ui_CoevolutionByPairs()
        self.__ui.setupUi(self)
        self.__result = result
        self.__msa = msa
        self.__visualization_structure : Optional[VisualizationStructure] = None
        self.__color_selector = pymol_residues_color_selector(
            self.__ui.resultsTable,
            self.__residue_selector
        )
        self.__distance_calculator = None
        self.__ui.resultsTable.itemSelectionChanged.connect(self.__update_details_table)

        # Add validators to the filter boxes, set them to inital values and
        # add the event handler for the filter button
        conf_validator = QDoubleValidator(0,1,2)
        self.__ui.confidenceMin.setValidator(conf_validator)
        self.__ui.confidenceMin.setText("0.5")
        self.__ui.confidenceMax.setValidator(conf_validator)
        self.__ui.confidenceMax.setText("1.0")

        dist_validator = QDoubleValidator(0, 1000, 2)
        self.__ui.distanceMin.setValidator(dist_validator)
        self.__ui.distanceMin.setText("0.0")
        self.__ui.confidenceMax.setValidator(dist_validator)
        self.__ui.distanceMax.setText("8.0")

        self.__ui.filterButton.clicked.connect(self.__render_results_table)

        # Call this after everything. Python is a dangerous language that
        # defines fields the moment they are assigned. Static type-checkers
        # will have difficulties copeing with such madness
        self.__update_results_table(result, msa)

    def __get_result_data(self, index : QModelIndex) -> ResultRoleType:
        return index.data(RESULT_ROLE)

    def __residue_selector(self, table : QTableWidget, index: QModelIndex) -> Optional[ResidueWithIntensity]:

        if self.__visualization_structure is None:
            return None

        result_data = self.__get_result_data(index)

        return ResidueWithIntensity(
            residues=set([
                result_data.pos_1_fallback,
                result_data.pos_2_fallback
            ]),
            intensity=result_data.entry.confidence
        )

    def set_visualization_structure(self, structure_sele: StructureSelection, link: List[int]):
        self.__visualization_structure = VisualizationStructure(
            structure=structure_sele,
            link = link,
            sequence_index=structure.get_pdb_sequence_index(structure_sele)
        )
        self.__color_selector.set_selection(structure_sele.selection)
        self.__distance_calculator = ResidueDistanceCalculator(structure_sele.selection)
        
        # This should be called last. That's what one gets when using imperative dynamic
        # languages like python.
        self.__render_results_table()

    @pyqtSlot()
    def __render_results_table(self):
        self.__update_results_table(self.__result, self.__msa)

    def __distance_in_visualization_structure(self, result: AcpsicovEntry) -> Optional[float]:

        if self.__distance_calculator is not None and self.__visualization_structure is not None:
            link = self.__visualization_structure.link
            distance = \
                self.__distance_calculator.distance_zero_index(
                    link[result.position1],
                    link[result.position2]
                )
            return round(distance, 2) if distance is not None else None
        else:
            return None

    def __set_result_data(self, item : QTableWidgetItem, entry : ResultRoleType):
        item.setData(RESULT_ROLE, entry)

    def __get_filtered_results(self, results: AcpsicovResult) -> Iterable[ResultItemModel]:

        min_conf = float(self.__ui.confidenceMin.text())
        max_conf = float(self.__ui.confidenceMax.text())
        min_dist = float(self.__ui.distanceMin.text())
        max_dist = float(self.__ui.distanceMax.text())

        for result in results.entries:

            if result.confidence < min_conf or result.confidence > max_conf:
                continue

            distance = self.__distance_in_visualization_structure(result)

            if distance is not None and (distance < min_dist or distance > max_dist):
                continue

            if self.__visualization_structure is not None:
                visualization_structure = self.__visualization_structure

                yield ResultItemModel(
                    pymol_position_1 = visualization_structure.msa_to_structure_pymol_posiiton(result.position1),
                    pymol_position_2 = visualization_structure.msa_to_structure_pymol_posiiton(result.position2),
                    entry = result,
                    distance=distance
                )
            else:
                yield ResultItemModel(
                    pymol_position_1=None,
                    pymol_position_2=None,
                    entry = result,
                    distance=distance
                )

    def __update_results_table(self, results : AcpsicovResult, msa : Msa):

        table = self.__ui.resultsTable
        table.clear()

        entries = list(self.__get_filtered_results(results))
        table.setRowCount(len(entries))
        
        visualization_structure = self.__visualization_structure

        if visualization_structure:
            structure_name = visualization_structure.structure.structure_name[0:20]
        else:
            structure_name = NA_VALUE

        columns = [
            "Position 1 (MSA)",
            f"Position 1 ({structure_name})",
            f"Residue 1 ({structure_name})",
            f"Position 2 (MSA)",
            f"Position 2 ({structure_name})",
            f"Residue 2 ({structure_name})",
            "Confidence",
            f"Distance ({structure_name}, Å)"
        ]
        table.setColumnCount(len(columns))
        table.setHorizontalHeaderLabels(columns)

        def get_at_position(msa_pos : int) -> Tuple[str, str]:
            if self.__visualization_structure is not None:
                return (
                    str(self.__visualization_structure.msa_to_structure_pymol_posiiton(msa_pos)),
                    self.__visualization_structure.msa_to_structure_residue(msa_pos)
                )
            else:
                return (NA_VALUE, NA_VALUE)

        for i,result in enumerate(entries):

            def __with_data__(text : str):

                item = QTableWidgetItem(text)
                self.__set_result_data(item, result)
                return item

            (pos1_structure, res1_structure) = get_at_position(result.entry.position1)
            (pos2_structure, res2_structure) = get_at_position(result.entry.position2)

            table.setItem(i, 0, __with_data__(str(result.entry.position1)))
            table.setItem(i, 1, __with_data__(pos1_structure))
            table.setItem(i, 2, __with_data__(res1_structure))
            table.setItem(i, 3, __with_data__(str(result.entry.position2)))
            table.setItem(i, 4, __with_data__(pos2_structure))
            table.setItem(i, 5, __with_data__(res2_structure))
            table.setItem(i, 6, __with_data__(str(result.entry.confidence)))
            table.setItem(i, 7, __with_data__(result.distance_str))

    @pyqtSlot()
    def __update_details_table(self):
        
        table = self.__ui.detailsTable
        table.clear()
        selection = self.__ui.resultsTable.selectedIndexes()

        if len(selection) != 1:
            return

        selected_result = selection[0]
        entry = self.__get_result_data(selected_result).entry
        pairs = enumerate_pairs(self.__msa, entry.position1 - 1, entry.position2 - 1)
        headers = list(sorted(pairs.keys(), key = lambda k: pairs[k], reverse=True))

        table.setColumnCount(len(headers))
        table.setRowCount(1)
        table.setHorizontalHeaderLabels(f"({p1}, {p2})" for p1, p2 in headers)

        for i,pair in enumerate(headers):
            table.setItem(0, i, QTableWidgetItem(str(pairs[pair])))