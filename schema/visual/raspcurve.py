from PyQt5.QtWidgets import QMainWindow, QSizePolicy, QTableWidgetItem, QTableWidget, QVBoxLayout, QWidget
import pymol

from .SchemaContext import SchemaContext
from ..models import SchemaResult

dummy_results = [(46.8400,140.1468,[112,180,229,452,521]) \
                ,(48.1600,140.7471,[112,180,230,448,521])] 

colors = ["c%i00" % x for x in range(1,9)] + ["c%i50" % x for x in range(1,9)]

def get_color(item, options = colors):
    return options[item % len(colors)]

class RaspCurveViewer(QWidget):

    def __init__(self, schema_context, *args, **kwargs):
        super(RaspCurveViewer, self).__init__(*args, **kwargs)
        self.__schema_context = schema_context
        layout = QVBoxLayout(self)
        self.__results_table = QTableWidget()
        layout.addWidget(self.__results_table)
        self.__results_table.cellActivated.connect(self.__on_cell_activated)

    def __reset(self):
        self.__results_table.clearContents()
        self.__results = None
        pass

    def __on_cell_activated(self, row, column):
        item = self.__results[row]

        for i in item[2]:
            s = item[i]
            e = item[i+1]
            pymol.cmd.color(get_color(i), "resi %i-%i" % (s, e))
        
    def set_results(self, results : SchemaResult):

        self.__reset()
        
        if len(results.results) < 1:
            self.__schema_context.raise_error_message("There must be at least one result.")
            return

        self.__results = results
        self.__results_table.setColumnCount(3)
        self.__results_table.setRowCount(len(results))

        for (row, (energy, mutations, points)) in enumerate(results.results):
            self.__results_table.setItem(row, 0, QTableWidgetItem(str(energy)))
            self.__results_table.setItem(row, 1, QTableWidgetItem(str(mutations)))
            self.__results_table.setItem(row, 2, QTableWidgetItem(str(points)))        
