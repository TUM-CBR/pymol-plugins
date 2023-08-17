from os import path
from PyQt5 import QtCore
from PyQt5.QtCore import QUrl, Qt, pyqtSlot
from PyQt5.QtGui import QColor, QDesktopServices
from PyQt5.QtWidgets import QListWidgetItem, QTableWidgetItem, QWidget
import pymol
from typing import Dict, List

from ...core.Qt.QtWidgets import get_qtable_content, open_copy_context_menu
from ...core.pymol import structure
from ...core.pymol.structure import StructureSelection
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

    @property
    def __result_selector(self) -> StructureSelection:
        return StructureSelection(
            structure_name=self.__result.structure_name,
            chain_name=None,
            segment_identifier=None
        )

    def __show_result_sequences(self, parents_msa : Dict[str, str]):

        table = self.__ui.alignmentViewer
        values = list(parents_msa.values())
        residues = len(values[0])

        table.clear()
        table.setRowCount(len(parents_msa))
        table.setColumnCount(residues + 1)
        table.setHorizontalHeaderLabels(
            ["Parent"] + [str(i) for i in range(1, residues + 1)]
        )

        for i in range(residues):
            for j, (k,v) in enumerate(parents_msa.items()):

                # Column 0 will contain the name of the parent
                cx = i + 1
                table.setItem(j, 0, QTableWidgetItem(k))

                residue_item = QTableWidgetItem(v[i])
                residue_item.setFlags(Qt.ItemIsEnabled)
                table.setItem(j, cx, residue_item)

        table.resizeColumnsToContents()

    def __set_result(self, result : SchemaResult):

        pymol.cmd.load(result.pdb)

        self.__result = result
        self.__result_items = result.load_results()

        results_viewer = self.__ui.resultsViewer
        results_viewer.clearContents()
        results_viewer.setRowCount(len(self.__result_items))
        mappings = structure.get_pdb_sequence_index(self.__result_selector)

        parents_msa = self.__result.parents_msa()
        self.__show_result_sequences(parents_msa)

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
            self.__ui.resultsViewer.setItem(row, 3, QTableWidgetItem(" ".join(map(str, item.shuffling_points(mappings)))))

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

    def __color_table_sequences(self, shuffling_points : List[int]):
        table = self.__ui.alignmentViewer

        previous_position = 1
        last_position = table.columnCount() - 1
        for i,position in enumerate(shuffling_points + [last_position]):
            for rx in range(previous_position, position + 1):
                for row in range(0, table.rowCount()):
                    target = table.item(row, rx)
                    assert target, "Bug in the code. We should always get a valid cell"
                    target.setBackground(QColor(*get_qt_color(i)))
            previous_position = position + 1

    def __select_item(self, row):
        item = self.__result_items[row]
        offset = structure.get_pdb_sequence_index(self.__result_selector)
        shuffling_points = item.shuffling_points(offset)

        # Ensure that we color the whole structure so the last fragment
        # will remain with this color
        pymol.cmd.color(
            get_color(len(shuffling_points)),
            "model %s" % self.__result.structure_name
        )

        e = 0

        for (i,loc) in enumerate(shuffling_points):
            s = e + 1
            e = loc
            sele = "(model %s) and (resi %i-%i)" % (self.__result.structure_name, s, e)
            pymol.cmd.color(get_color(i), sele)

        self.__color_table_sequences(item.msa_shuffling_points)

distinct_colors = [
    (255, 0, 0),
    (0, 255, 0),
    (0, 0, 255),
    (255, 255, 0),
    (0, 255, 255),
    (255, 0, 255),
    (128, 128, 0),
    (0, 128, 128),
    (128, 0, 128),
    (128, 0, 0),
    (0, 128, 0),
    (0, 0, 128),
    (192, 192, 192),
    (128, 128, 128),
    (255, 255, 255)
]

def to_hex(i : int, desired_length = 2) -> str:
    value = hex(i).replace("0x", "")
    padding = desired_length - len(value)
    return "0"*padding + value

colors = [
    "0x" + "".join(map(to_hex, [r,g,b]))
    for r,g,b in distinct_colors
]

def get_qt_color(item : int, options = distinct_colors) -> List[int]:
    return list(options[item % len(colors)]) + [100]

def get_color(item, options = colors):
    return options[item % len(colors)]