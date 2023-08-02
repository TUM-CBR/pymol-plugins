from os import path
from PyQt5 import QtCore
from PyQt5.QtCore import QUrl, Qt, pyqtSlot
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtWidgets import QListWidgetItem, QTableWidgetItem, QWidget
import pymol

from ...core.Qt.QtWidgets import get_qtable_content, open_copy_context_menu
from ...core.pymol import structure
from ...clustal import msa

from ..SchemaContext import SchemaContext
from ..SchemaResult import SchemaResult
from ..SchemaTaskManager import SchemaTaskManager

from .Ui_SchemaSelectorWidget import Ui_SchemaSelectorWidget

class SchemaSelectorWidget(QWidget):

    RESULT_DATA_ROLE = 1
    MAIN_COLUMN_HEADERS = ["<E>", "<M>", "Positions (MSA)", "Positions (Structure)"]

    def __init__(self, schema_context, manager : SchemaTaskManager, *args, **kwargs):
        super(SchemaSelectorWidget, self).__init__(*args, **kwargs)
        self.__schema_context : SchemaContext = schema_context
        self.__ui = Ui_SchemaSelectorWidget()
        self.__ui.setupUi(self)
        self.__manager = manager
        self.__manager.subscribe_results_updated(self.__on_results_updated)
        self.__ui.resultsViewer.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.__ui.resultsViewer.customContextMenuRequested.connect(self.__on_results_context_menu)
        self.__ui.resultsList.itemDoubleClicked.connect(self.__on_item_double_clicked)

    @pyqtSlot(QtCore.QPoint)
    def __on_results_context_menu(self, pos):
        open_copy_context_menu(self.__ui.resultsViewer, pos)

    def __on_results_updated(self, results):
        self.__ui.resultsList.clear()

        for result in results:
            item = QListWidgetItem(result.name, self.__ui.resultsList)
            item.setData(SchemaSelectorWidget.RESULT_DATA_ROLE, result)

    def __set_result(self, result : SchemaResult):

        pymol.cmd.load(result.pdb)

        self.__result = result
        self.__result_items = result.load_results()

        results_viewer = self.__ui.resultsViewer
        results_viewer.clearContents()
        results_viewer.setRowCount(len(self.__result_items))
        offset = structure.get_structure_offset(self.__result.structure_name)

        parents_msa = self.__result.parents_msa()

        # We have the four fixed columns for <E> <M> "MSA positions" and "Structure Positions"
        # We add one column per sequence in the results
        columns = SchemaSelectorWidget.MAIN_COLUMN_HEADERS + ["Positions (%s)" % k for k in parents_msa.keys()]
        results_viewer.setColumnCount(len(columns))
        results_viewer.setHorizontalHeaderLabels(columns)

        def get_sequence_position(i : int, seq_name : str) -> int:
            return len(msa.clean_msa_blanks(parents_msa[seq_name][0:i]))

        for (row, item) in enumerate(self.__result_items):
            self.__ui.resultsViewer.setItem(row, 0, QTableWidgetItem(item.energy))
            self.__ui.resultsViewer.setItem(row, 1, QTableWidgetItem(item.mutations))
            self.__ui.resultsViewer.setItem(row, 2, QTableWidgetItem(" ".join(map(str, item.msa_shuffling_points))))
            self.__ui.resultsViewer.setItem(row, 3, QTableWidgetItem(" ".join(map(str, item.shuffling_points(offset)))))

            for (i, seq_name) in enumerate(parents_msa.keys()):
                ix = i + 4
                positions = [str(get_sequence_position(sp, seq_name)) for sp in item.msa_shuffling_points]
                results_viewer.setItem(row, ix, QTableWidgetItem(" ".join(positions)))

        results_viewer.resizeColumnsToContents()

        directory = path.dirname(result.results_file)
        results_with_positions = path.join(directory, "results_with_pdb_positions.txt")

        if not path.exists(results_with_positions):
            with open(results_with_positions, 'w') as results_stream:
                results_stream.write(get_qtable_content(self.__ui.resultsViewer, lambda _1,_2: True))

    @pyqtSlot(QListWidgetItem)
    def __on_item_double_clicked(self, item : QListWidgetItem):
        result : SchemaResult = item.data(SchemaSelectorWidget.RESULT_DATA_ROLE)
        QDesktopServices.openUrl(QUrl.fromLocalFile(path.dirname(result.results_file)))

    @pyqtSlot(QListWidgetItem)
    def on_resultsList_itemClicked(self, item : QListWidgetItem):
        self.__set_result(item.data(SchemaSelectorWidget.RESULT_DATA_ROLE))

    def on_resultsViewer_cellActivated(self, row, column):
        self.__select_item(row)

    def on_resultsViewer_cellClicked(self, row, column):
        self.__select_item(row)

    def __select_item(self, row):
        item = self.__result_items[row]
        offset = structure.get_structure_offset(self.__result.structure_name)
        shuffling_points = item.shuffling_points(offset)

        # Ensure that we color the whole structure so the last fragment
        # will remain with this color
        pymol.cmd.color(
            get_color(len(shuffling_points)),
            "model %s" % self.__result.structure_name
        )

        e = 0

        for (i,loc) in enumerate(shuffling_points):
            s = e
            e = loc
            sele = "(model %s) and (resi %i-%i)" % (self.__result.structure_name, s, e)
            pymol.cmd.color(get_color(i), sele)

colors = ["c%i00" % x for x in range(1,9)] + ["c%i50" % x for x in range(1,9)]

def get_color(item, options = colors):
    return options[item % len(colors)]