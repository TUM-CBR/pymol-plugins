from PyQt5.QtWidgets import QMainWindow

from .SchemaContext import SchemaContext
from .visual.SchemaInstanceWidget import SchemaInstanceWidget

schema_ctx = SchemaContext()
mw = QMainWindow()
schema_ctx = SchemaContext()
tw = SchemaInstanceWidget(schema_ctx)
mw.setCentralWidget(tw)
mw.show()

#tw.set_results(dummy_results)
