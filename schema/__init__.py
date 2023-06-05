from PyQt5.QtWidgets import QMainWindow, QPushButton, QTableWidget

from .visual.SchemaContext import SchemaContext
from .visual.raspcurve import dummy_results, RaspCurveViewer
from .visual.SchemaRunner import SchemaRunner

mw = QMainWindow()
schema_ctx = SchemaContext()
# tw = QPushButton(text="hello")
tw = RaspCurveViewer(schema_ctx)
tw = SchemaRunner(schema_ctx)
mw.setCentralWidget(tw)
mw.show()

#tw.set_results(dummy_results)
