from typing import List

from PyQt5.QtWidgets import QTableWidgetItem

from ...acpsicov.main import AcpsicovResult
from ...core.pymol.structure import StructureSelection
from ...support.msa import Msa

from .CoevolutionResultViewerBase import CoevolutionResultViewerBase
from .Ui_CoevolutionByPairs import Ui_CoevolutionByPairs

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

        self.__update_results_table(result, msa)

    def set_visualization_structure(self, structure: StructureSelection,link: List[int]):
        return

    def __update_results_table(self, results : AcpsicovResult, msa : Msa):

        table = self.__ui.resultsTable
        table.clear()

        table.setRowCount(len(results.entries))
        table.setColumnCount(3)

        table.setHorizontalHeaderLabels([
            "Position 1 (MSA)",
            "Position 2 (MSA)",
            "Confidence"
        ])


        for i,result in enumerate(results.entries):
            table.setItem(i, 0, QTableWidgetItem(str(result.position1)))
            table.setItem(i, 1, QTableWidgetItem(str(result.position2)))
            table.setItem(i, 2, QTableWidgetItem(str(result.confidence)))