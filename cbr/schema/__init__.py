def __init_plugin__(app=None):

    from PyQt5.QtWidgets import QMainWindow

    from .SchemaContext import SchemaContext
    from .visual.SchemaInstanceWidget import SchemaInstanceWidget

    schema_ctx = SchemaContext()
    mw = QMainWindow()
    schema_ctx = SchemaContext()
    tw = SchemaInstanceWidget(schema_ctx)
    mw.setCentralWidget(tw)

    from pymol.plugins import addmenuitemqt
    addmenuitemqt('CBR Tools', mw.show)
