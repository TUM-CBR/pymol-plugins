from PyQt5.QtWidgets import QListWidgetItem, QMainWindow, QSizePolicy, QTableWidgetItem, QTableWidget, QVBoxLayout, QWidget
import pymol

from ..SchemaResult import SchemaResult

from .SchemaContext import SchemaContext
from .Ui_SchemaSelectorWidget import Ui_SchemaSelectorWidget
from ..SchemaTaskManager import SchemaTaskManager

class SchemaSelectorWidget(QWidget):

    RESULT_DATA_ROLE = 1

    def __init__(self, schema_context, manager : SchemaTaskManager, *args, **kwargs):
        super(SchemaSelectorWidget, self).__init__(*args, **kwargs)
        self.__schema_context = schema_context
        self.__ui = Ui_SchemaSelectorWidget()
        self.__ui.setupUi(self)
        self.__manager = manager
        self.__manager.subscribe_results_updated(self.__on_results_updated)

    def __on_results_updated(self, results):
        self.__ui.resultsList.clear()

        for result in results:
            item = QListWidgetItem(result.name, self.__ui.resultsList)
            item.setData(SchemaSelectorWidget.RESULT_DATA_ROLE, result)

    def __set_result(self, result : SchemaResult):

        self.__result = result
        self.__result_items = result.load_results()

        self.__ui.resultsViewer.clearContents()
        self.__ui.resultsViewer.setColumnCount(3)
        self.__ui.resultsViewer.setRowCount(len(self.__result_items))

        for (row, item) in enumerate(self.__result_items):
            self.__ui.resultsViewer.setItem(row, 0, QTableWidgetItem(item.energy))
            self.__ui.resultsViewer.setItem(row, 1, QTableWidgetItem(item.mutations))
            self.__ui.resultsViewer.setItem(row, 2, QTableWidgetItem(str(item.shuffling_points)))

        pymol.cmd.load(result.pdb)

    def on_resultsList_itemClicked(self, item : QListWidgetItem):
        self.__set_result(item.data(SchemaSelectorWidget.RESULT_DATA_ROLE))

    def on_resultsViewer_cellActivated(self, row, column):
        self.__select_item(row)

    def on_resultsViewer_cellClicked(self, row, column):
        self.__select_item(row)

    def __select_item(self, row):
        item = self.__result_items[row]

        e = 0
        for (i,loc) in enumerate(item.shuffling_points):
            s = e
            e = loc
            sele = "(model %s) and (resi %i-%i)" % (self.__result.structure_name, s, e)
            pymol.cmd.color(get_color(i), sele)

dummy_results = [(46.8400,140.1468,[112,180,229,452,521]) \
                ,(48.1600,140.7471,[112,180,230,448,521])] 

colors = ["c%i00" % x for x in range(1,9)] + ["c%i50" % x for x in range(1,9)]

def get_color(item, options = colors):
    return options[item % len(colors)]