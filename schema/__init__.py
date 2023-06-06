from PyQt5.QtWidgets import QMainWindow, QPushButton, QTableWidget

from .visual.SchemaContext import SchemaContext
from .visual.SchemaInstanceWidget import SchemaInstanceWidget

mw = QMainWindow()
schema_ctx = SchemaContext()
# tw = QPushButton(text="hello")
tw = SchemaInstanceWidget(schema_ctx)
mw.setCentralWidget(tw)
mw.show()

#tw.set_results(dummy_results)
