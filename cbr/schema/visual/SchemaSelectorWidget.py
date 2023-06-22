from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QAction, QApplication, QMenu, QListWidgetItem, QTableWidgetItem, QWidget
import pymol

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
        self.__schema_context = schema_context
        self.__ui = Ui_SchemaSelectorWidget()
        self.__ui.setupUi(self)
        self.__manager = manager
        self.__manager.subscribe_results_updated(self.__on_results_updated)
        self.__ui.resultsViewer.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.__ui.resultsViewer.customContextMenuRequested.connect(self.__on_results_context_menu)

    def __on_results_context_menu(self, pos):
        context_menu = QMenu(self)
        copy_action = QAction("Copy", self)
        copy_action.triggered.connect(self.__copyTable)
        context_menu.addAction(copy_action)
        context_menu.exec_(self.__ui.resultsViewer.mapToGlobal(pos))

    def __copyTable(self):
        table_widget = self.__ui.resultsViewer
        selected_indexes = table_widget.selectedIndexes()
        if len(selected_indexes) > 0:
            selected_cells = set()
            for index in selected_indexes:
                selected_cells.add((index.row(), index.column()))

            copied_data = ""
            for row in range(table_widget.rowCount()):
                for column in range(table_widget.columnCount()):
                    if (row, column) in selected_cells:
                        item = table_widget.item(row, column)
                        copied_data += str(item.text()) + "\t"
                copied_data += "\n"

            clipboard = QApplication.clipboard()
            clipboard.setText(copied_data)

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

    def on_resultsList_itemClicked(self, item : QListWidgetItem):
        self.__set_result(item.data(SchemaSelectorWidget.RESULT_DATA_ROLE))

    def on_resultsViewer_cellActivated(self, row, column):
        self.__select_item(row)

    def on_resultsViewer_cellClicked(self, row, column):
        self.__select_item(row)

    def __select_item(self, row):
        item = self.__result_items[row]
        offset = self.__get_offset(self.__result.structure_name)
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