import os
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt, pyqtSlot
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QVBoxLayout, QWidget
from typing import Optional

from ...core.Qt.QtCore import run_in_thread
from ...core.Qt.QtWidgets import progress_manager

from ..data import DesignPrimersResults, PrimerResult
from .. import operations

from .PrimerResultViewer import PrimerResultContext, PrimerResultViewer
from .Ui_PrimerViewer import Ui_PrimerViewer

K_ITEM_DATA_ROLE = 1

def render_tm(tm : float):
    return f"{round(tm, ndigits=1)}° C"

class PrimersDataModel(QAbstractTableModel):

    def __init__(self, primers : DesignPrimersResults):
        super().__init__()
        self.__primers = primers
        self.__columns = [
            "Position",
            "New Residue",
            "Tm (Sequence)",
            "Tm (Left Primer)",
            "Tm (Right Primer)",
            "Codon",
            "Left Primer",
            "Right Primer"
        ]

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            headers = self.__columns
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def rowCount(self, parent=None):
        return len(self.__primers.primers)

    def columnCount(self, parent=None):
        return len(self.__columns)
    
    def data(self, index: QModelIndex, role=Qt.DisplayRole):

        if not index.isValid():
            return None

        primer = self.__primers.primers[index.row()]

        if role == Qt.DisplayRole:
            cols = [
                str(primer.position),
                str(primer.amino_acid),
                render_tm(primer.tm_all),
                render_tm(primer.tm_left),
                render_tm(primer.tm_right),
                str(primer.inner_seq),
                str(primer.left_primer),
                str(primer.right_primer)
            ]
            self.setItemData(index, {K_ITEM_DATA_ROLE: primer})
            return cols[index.column()]
        elif role == K_ITEM_DATA_ROLE:
            return primer
        else:
            return None


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

        # Populate the tables with some initial values
        self.__on_select_primers()

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
        return self.__ui.primersTable.model().data(index, K_ITEM_DATA_ROLE)

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
                source = scoped_results
            )
        )

    @pyqtSlot(object)
    def __render_results(self, results : DesignPrimersResults):
        self.__render_primers_table(results)
        self.__primer_viewer.set_result(
            PrimerResultContext(
                result = None,
                source = results
            )
        )

    def __render_primers_table(self, results : DesignPrimersResults):

        self.__scoped_results = results
        table = self.__ui.primersTable
        model = PrimersDataModel(results)
        table.setModel(model)
        table.selectionModel().selectionChanged.connect(self.__on_item_selection_changed)
        table.resizeColumnsToContents()