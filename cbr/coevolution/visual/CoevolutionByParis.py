from PyQt5.QtCore import QModelIndex
from typing import List, NamedTuple, Optional, Tuple

from PyQt5.QtWidgets import QTableWidget, QTableWidgetItem

from ...acpsicov.main import AcpsicovEntry, AcpsicovResult
from ...core.pymol.structure import StructureSelection
from ...core.pymol import structure
from ...core.pymol.visual.PymolResidueResultsTable import pymol_residues_color_selector, ResidueWithIntensity
from ...support.msa import Msa

from .CoevolutionResultViewerBase import CoevolutionResultViewerBase
from .Ui_CoevolutionByPairs import Ui_CoevolutionByPairs

class VisualizationStructure(NamedTuple):
    structure : StructureSelection
    link : List[int]

RESULT_ROLE = 1

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

        self.__update_results_table(result, msa)

    def __residue_selector(self, table : QTableWidget, index: QModelIndex) -> Optional[ResidueWithIntensity]:

        if self.__visualization_structure is None:
            return None


        result_data : AcpsicovEntry = index.data(RESULT_ROLE)

        return ResidueWithIntensity(
            residues=set([
                result_data.position1,
                result_data.position2
            ]),
            intensity=result_data.confidence
        )

    def set_visualization_structure(self, structure: StructureSelection,link: List[int]):
        self.__visualization_structure = VisualizationStructure(
            structure=structure,
            link = link
        )
        self.__color_selector.set_selection(structure.selection)
        self.__update_results_table(self.__result, self.__msa)

    def __update_results_table(self, results : AcpsicovResult, msa : Msa):

        table = self.__ui.resultsTable
        table.clear()

        table.setRowCount(len(results.entries))
        table.setColumnCount(7)
        
        visualization_structure = self.__visualization_structure

        if visualization_structure:
            structure_name = f" ({visualization_structure.structure.structure_name[0:20]})"
            structure_sequence_index = structure.get_pdb_sequence_index(visualization_structure.structure)
            structure_sequence_offset = min(structure_sequence_index.keys())
            structure_sequence = "".join(structure_sequence_index.values())
            link = visualization_structure.link
        else:
            structure_name = " (N/A)"
            structure_sequence = None
            structure_sequence_offset = 0
            link = None

        table.setHorizontalHeaderLabels([
            "Position 1 (MSA)",
            "Position 1" + structure_name,
            "Residue 1" + structure_name,
            "Position 2 (MSA)",
            "Position 2" + structure_name,
            "Residue 2" + structure_name,
            "Confidence",
        ])

        def get_at_position(msa_pos : int) -> Tuple[str, str]:
            if link is not None and structure_sequence is not None:
                structure_pos = link[msa_pos]
                return (
                    str(structure_pos + structure_sequence_offset),
                    structure_sequence[structure_pos]
                )
            else:
                return ("N/A", "N/A")

        for i,result in enumerate(results.entries):

            def __with_data__(text : str):

                if link is not None:
                    position1 = link[result.position1] + structure_sequence_offset
                    position2 = link[result.position2] + structure_sequence_offset
                else:
                    position1 = result.position1
                    position2 = result.position2

                item = QTableWidgetItem(text)
                item.setData(
                    RESULT_ROLE,
                    result._replace(
                        position1 = position1,
                        position2 = position2
                    )
                )
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