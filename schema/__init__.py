from PyQt5.QtWidgets import QMainWindow

from .visual.SchemaContext import SchemaContext
from .visual.SchemaInstanceWidget import SchemaInstanceWidget

if not SchemaContext.is_loaded():
    mw = QMainWindow()
    schema_ctx = SchemaContext()
    # tw = QPushButton(text="hello")
    tw = SchemaInstanceWidget(schema_ctx)
    mw.setCentralWidget(tw)
    mw.show()

#tw.set_results(dummy_results)
