from PyQt5.QtWidgets import QMainWindow, QSizePolicy, QTableWidgetItem, QTableWidget, QVBoxLayout, QWidget
import pymol

dummy_results = [(46.8400,140.1468,112,180,229,452,521) \
                ,(48.1600,140.7471,112,180,230,448,521)] 

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
        # self.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        # self.__results_table.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        self.__results_table.cellActivated.connect(self.__on_cell_activated)

    def __reset(self):
        self.__results_table.clearContents()
        self.__results = None
        pass

    def __on_cell_activated(self, row, column):
        print("le click my boy")
        item = self.__results[row]

        for i in range(0, len(item) - 1):
            if i > 1:
                s = item[i]
                e = item[i+1]
                pymol.cmd.color(get_color(i), "resi %i-%i" % (s, e))
        
    def set_results(self, results):

        self.__reset()
        
        if len(results) < 1 \
           and len(results[0]) > 3 \
           and any([len(r) != len(results[0]) for r in results]):
            self.__schema_context.raise_error_message("There must be at least one result and all results must have the same number of shuffle points")
            return

        self.__results = results
        self.__results_table.setColumnCount(len(results[0]))
        self.__results_table.setRowCount(len(results))

        for (row, result) in enumerate(results):
            for (column, value) in enumerate(result):
                self.__results_table.setItem(row, column, QTableWidgetItem(str(value)))
