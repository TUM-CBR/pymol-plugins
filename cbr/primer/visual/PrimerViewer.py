import os
from PyQt5.QtCore import QModelIndex, pyqtSlot
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QTableWidgetItem, QVBoxLayout, QWidget
from typing import Dict, Optional

from ...core.Qt.QtCore import run_in_thread
from ...core.Qt.QtWidgets import progress_manager

from ..data import DesignPrimersResults, PrimerResult
from .. import operations

from .PrimerResultViewer import PrimerResultContext, PrimerResultViewer
from .Ui_PrimerViewer import Ui_PrimerViewer

K_ITEM_DATA_ROLE = 1

DesignPrimersFilteredResult = Dict[str, PrimerResult]

class PrimerViewer(QWidget):

    def __init__(
        self,
        database : str,
        delete_when_closed = False
    ) -> None:
        super().__init__()

        self.__ui = Ui_PrimerViewer()
        self.__ui.setupUi(self)
        self.__database = database
        self.__scoped_results : Optional[DesignPrimersResults] = None
        self.__delete_when_closed = delete_when_closed

        layout = QVBoxLayout()
        self.__primer_viewer = PrimerResultViewer(self.__ui.resultsWidget)
        layout.addWidget(self.__primer_viewer)
        self.__ui.resultsWidget.setLayout(layout)

        self.__ui.tmEdit.setValidator(QDoubleValidator(0, 100, 2))
        self.__ui.tmEdit.setText("65.00")
        self.__ui.tmWeightEdit.setValidator(QDoubleValidator(0,100,2))
        self.__ui.tmWeightEdit.setText("1.00")
        self.__ui.tmPrimersEdit.setValidator(QDoubleValidator(0,100,2))
        self.__ui.tmPrimersEdit.setText("45.00")
        self.__ui.tmPrimersWeightEdit.setValidator(QDoubleValidator(0,100,2))
        self.__ui.tmPrimersWeightEdit.setText("1.00")
        self.__ui.tmDeltaWeightEdit.setValidator(QDoubleValidator(0,100,2))
        self.__ui.tmDeltaWeightEdit.setText("1.00")
        self.__ui.selectButton.clicked.connect(self.__on_select_primers)

        self.__progress = progress_manager(
            self.__ui.filterProgress,
            self.__ui.selectButton
        )
        self.__progress.with_default_error_handler(self)
        self.__progress.on_result.connect(self.__render_results)
        self.__ui.primersTable.itemSelectionChanged.connect(self.__on_item_selection_changed)

    def __del__(self):

        if not self.__delete_when_closed:
            return

        try:
            os.remove(self.__database)
        except Exception:
            # Ignore if file could not be deleted
            pass

    @pyqtSlot()
    def __on_select_primers(self):

        self.__progress.watch_progress(
            self.__select_primers(
                float(self.__ui.tmEdit.text()),
                float(self.__ui.tmWeightEdit.text()),
                float(self.__ui.tmPrimersEdit.text()),
                float(self.__ui.tmPrimersWeightEdit.text()),
                float(self.__ui.tmDeltaWeightEdit.text())
            )
        )

    @run_in_thread
    def __select_primers(
        self,
        tm : float,
        wTm : float,
        pTm : float,
        wPTm : float,
        wTmDelta : float
    ) -> DesignPrimersResults:
        return operations.query_best_primers(
            self.__database,
            tm,
            wTm,
            pTm,
            wPTm,
            wTmDelta
        )

    def __get_primer_from_index(self, index : QModelIndex) -> PrimerResult:
        return index.data(K_ITEM_DATA_ROLE)

    @pyqtSlot()
    def __on_item_selection_changed(self):
        table = self.__ui.primersTable
        indexes = table.selectedIndexes()

        if len(indexes) <= 0:
            return
        ix = min(indexes, key = lambda item: item.row())
        item = self.__get_primer_from_index(ix)

        scoped_results = self.__scoped_results
        assert scoped_results, "Primer selected withoth a plasmid in socpe! This is a bug."

        self.__primer_viewer.set_result(
            PrimerResultContext(
                result = item,
                source=scoped_results
            )
        )

    @pyqtSlot(object)
    def __render_results(self, results : DesignPrimersResults):

        self.__scoped_results = results

        table = self.__ui.primersTable
        total_rows = len(results.primers)
        table.clearContents()
        table.setRowCount(total_rows)
        table.setColumnCount(7)
        columns = [
            "Position",
            "New Residue",
            "Tm (Sequence)",
            "Tm (Left Primer)",
            "Tm (Right Primer)",
            "Codon",
            "Left Primer",
            "Right Primer"
        ]
        table.setHorizontalHeaderLabels(columns)

        for i,primer in enumerate(results.primers):

            def new_cell(args):
                item = QTableWidgetItem(args)
                item.setData(K_ITEM_DATA_ROLE, primer)
                return item

            cells = [
                new_cell(str(primer.position)),
                new_cell(str(primer.amino_acid)),
                new_cell(str(primer.tm_all)),
                new_cell(str(primer.tm_left)),
                new_cell(str(primer.tm_right)),
                new_cell(str(primer.inner_seq)),
                new_cell(str(primer.left_primer)),
                new_cell(str(primer.right_primer))
            ]

            for j,cell in enumerate(cells):
                table.setItem(i, j, cell)

        table.resizeColumnsToContents()