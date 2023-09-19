import os
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QFileDialog, QWidget
from typing import Dict, NamedTuple

from ...acpsicov import main as acpsicov

from ...clustal import Clustal
from ...coevolution.visual.CoevolutionResultViewerBase import CoevolutionResultViewerBase
from ...core.Context import Context
from ...support.msa import Msa
from ...support.msa.visual import as_alignment_link_selector

from .CoevolutionByParis import CoevolutionByPairs
from .Ui_CoevolutionViewer import Ui_CoevolutionViewer

class CoevolutionResult(NamedTuple):
    msa : Msa
    acpsicov_result : acpsicov.AcpsicovResult

class CoevolultionViewer(QWidget):

    def __init__(
        self,
        result : CoevolutionResult,
        context : Context
    ):
        super(CoevolultionViewer, self).__init__()

        self.__ui = Ui_CoevolutionViewer()
        self.__ui.setupUi(self)
        self.__result = result

        self.__alignment_link_selector = as_alignment_link_selector(
            self.__ui.structureCombo,
            self.__ui.refreshButton,
            self.__ui.sequenceCombo,
            result.msa,
            Clustal.get_clustal_from_context(context),
            self.__ui.selectionErrorLabel
        )

        self.__ui.exportResultsButton.clicked.connect(self.__on_export_results_clicked)

        self.__alignment_link_selector.struture_selected.connect(self.__on_structure_selection)

        self.__viewers : Dict[str, CoevolutionResultViewerBase] = {
            "View by Pairs" : CoevolutionByPairs(result.acpsicov_result, result.msa)
        }

        for name, widget in self.__viewers.items():
            self.__ui.tabWidget.addTab(
                widget,
                name
            )

    @pyqtSlot()
    def __on_export_results_clicked(self):

        file_name, _ = QFileDialog.getSaveFileName(
            self,
            "Save File",
            "",
            "ACPSICOV files (*.acpsicov)"
        )

        if os.path.splitext(file_name)[1].lower() != ".acpsicov":
            file_name = file_name + ".acpsicov"

        with open(file_name, 'w') as result_file:
            self.__result.acpsicov_result.write(result_file)

    @pyqtSlot()
    def __on_structure_selection(self):
        structure = self.__alignment_link_selector.structure
        link = self.__alignment_link_selector.relative_positions

        if structure is None or link is None:
            return

        for viewer in self.__viewers.values():
            viewer.set_visualization_structure(
                structure,
                link
            )
