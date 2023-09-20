from PyQt5.QtCore import QModelIndex, pyqtSlot
from typing import Dict, List, NamedTuple, Optional, Tuple

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
        return self.link[msa_pos]

    def msa_to_structure_pymol_posiiton(self, msa_pos: int) -> int:
        return self.msa_to_structure_position(msa_pos) + self.sequence_offset

    def msa_to_structure_residue(self, msa_pos: int) -> str:
        return self.sequence_index[self.msa_to_structure_pymol_posiiton(msa_pos)]

class ResultItemModel(NamedTuple):
    pymol_position_1 : Optional[int]
    pymol_position_2 : Optional[int]
    entry : AcpsicovEntry

    @property
    def pos_1_fallback(self) -> int:
        return self.pymol_position_1 if self.pymol_position_1 is not None else self.entry.position1

    @property
    def pos_2_fallback(self) -> int:
        return self.pymol_position_2 if self.pymol_position_2 is not None else self.entry.position2

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
        self.__update_results_table(self.__result, self.__msa)

    def __distance_in_visualization_structure(self, result: AcpsicovEntry):

        if self.__distance_calculator is not None and self.__visualization_structure is not None:
            link = self.__visualization_structure.link
            distance = \
                self.__distance_calculator.distance_zero_index(
                    link[result.position1],
                    link[result.position2]
                )
            return str(round(distance, 2)) if distance is not None else NA_VALUE
        else:
            return NA_VALUE

    def __set_result_data(self, item : QTableWidgetItem, entry : AcpsicovEntry):

        if self.__visualization_structure is not None:
            visualization_structure = self.__visualization_structure

            data : ResultRoleType = ResultItemModel(
                pymol_position_1 = visualization_structure.msa_to_structure_pymol_posiiton(entry.position1),
                pymol_position_2 = visualization_structure.msa_to_structure_pymol_posiiton(entry.position2),
                entry = entry
            )
        else:
            data : ResultRoleType = ResultItemModel(
                pymol_position_1=None,
                pymol_position_2=None,
                entry = entry
            )

        item.setData(RESULT_ROLE, data)

    def __update_results_table(self, results : AcpsicovResult, msa : Msa):

        table = self.__ui.resultsTable
        table.clear()

        table.setRowCount(len(results.entries))
        
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

        for i,result in enumerate(results.entries):

            def __with_data__(text : str):

                item = QTableWidgetItem(text)
                self.__set_result_data(item, result)
                return item

            (pos1_structure, res1_structure) = get_at_position(result.position1)
            (pos2_structure, res2_structure) = get_at_position(result.position2)

            table.setItem(i, 0, __with_data__(str(result.position1)))
            table.setItem(i, 1, __with_data__(pos1_structure))
            table.setItem(i, 2, __with_data__(res1_structure))
            table.setItem(i, 3, __with_data__(str(result.position2)))
            table.setItem(i, 4, __with_data__(pos2_structure))
            table.setItem(i, 5, __with_data__(res2_structure))
            table.setItem(i, 6, __with_data__(str(result.confidence)))
            table.setItem(i, 7, __with_data__(self.__distance_in_visualization_structure(result)))

    @pyqtSlot()
    def __update_details_table(self):
        
        table = self.__ui.detailsTable
        table.clear()
        selection = self.__ui.resultsTable.selectedIndexes()

        if len(selection) != 1:
            return

        selected_result = selection[0]
        entry = self.__get_result_data(selected_result).entry
        pairs = enumerate_pairs(self.__msa, entry.position1, entry.position2)

        headers = list(pairs.keys())
        table.setColumnCount(len(headers))
        table.setRowCount(1)
        table.setHorizontalHeaderLabels(f"({p1}, {p2})" for p1, p2 in headers)

        for i,pair in enumerate(headers):
            table.setItem(0, i, QTableWidgetItem(str(pairs[pair])))