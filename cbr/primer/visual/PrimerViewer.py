from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QColor, QDoubleValidator
from PyQt5.QtWidgets import QTableWidgetItem, QVBoxLayout, QWidget
from typing import Dict, Iterable, List

from ...core.Qt.QtCore import run_in_thread
from ...core.Qt.QtWidgets import progress_manager
from ...core.typing import from_just

from ..Operations import DesignPrimersResult, PrimerResult

from .PrimerResultViewer import PrimerResultContext, PrimerResultViewer
from .Ui_PrimerViewer import Ui_PrimerViewer

K_ITEM_DATA_ROLE = 1

DesignPrimersFilteredResult = Dict[str, PrimerResult]

class PrimerViewer(QWidget):

    def __init__(
        self,
        results : List[DesignPrimersResult],
        sequence : str
    ) -> None:
        super().__init__()

        self.__ui = Ui_PrimerViewer()
        self.__ui.setupUi(self)
        self.__sequence = sequence

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

        self.__results = results

        self.__progress = progress_manager(
            self.__ui.filterProgress,
            self.__ui.selectButton
        )
        self.__progress.with_default_error_handler(self)
        self.__progress.on_result.connect(self.__render_results)
        self.__ui.primersTable.itemSelectionChanged.connect(self.__on_item_selection_changed)

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
    ) -> List[DesignPrimersFilteredResult]:

        def select_best(primers : List[PrimerResult]) -> PrimerResult:

            primer = primers[0]

            for candidate in primers[1:]:
                primer = PrimerResult.choose_best(
                    primer,
                    candidate,
                    tm,
                    wTm,
                    pTm,
                    wPTm,
                    wTmDelta
                )
            
            return primer

        def select_primers_iterator() -> Iterable[DesignPrimersFilteredResult]:

            for result in self.__results:
                yield dict(
                    (resi, select_best(primers))
                    for resi, primers in result.items()
                )

        return list(select_primers_iterator())

    @pyqtSlot()
    def __on_item_selection_changed(self):
        table = self.__ui.primersTable
        indexes = table.selectedIndexes()

        if len(indexes) <= 0:
            return
        ix = min(indexes, key = lambda item: item.row())
        item = from_just(table.item(ix.row(), 0)).data(K_ITEM_DATA_ROLE)
        self.__primer_viewer.set_result(
            PrimerResultContext(
                result = item,
                original_sequence = self.__sequence
            )
        )

    @pyqtSlot(object)
    def __render_results(self, results : List[DesignPrimersFilteredResult]):
        table = self.__ui.primersTable
        total_rows = sum(
            len(result) for result in results
        )
        table.clearContents()
        table.setRowCount(total_rows)
        table.setColumnCount(7)
        columns = ["New Residue", "Tm (Sequence)", "Tm (Left Primer)", "Tm (Right Primer)", "Codon", "Left Primer", "Right Primer"]
        table.setHorizontalHeaderLabels(columns)

        i=0
        for result in results:
            for aa, primer in result.items():

                def new_cell(args):
                    item = QTableWidgetItem(args)
                    item.setData(K_ITEM_DATA_ROLE, primer)
                    return item

                table.setItem(i, 0, new_cell(aa))

                tm_all = primer.tm_all

                if tm_all:
                    tm_item = new_cell(str(primer.tm_all))
                else:
                    tm_item = new_cell(primer.tm_error or "")
                    tm_item.setBackground(QColor(100, 0, 0, 100))

                table.setItem(i, 1, tm_item)
                table.setItem(i, 2, new_cell(str(primer.tm_left)))
                table.setItem(i, 3, new_cell(str(primer.tm_right)))
                table.setItem(i, 4, new_cell(str(primer.inner_seq)))
                table.setItem(i, 5, new_cell(str(primer.left_primer)))
                table.setItem(i, 6, new_cell(str(primer.right_primer)))
                i+=1

        table.resizeColumnsToContents()