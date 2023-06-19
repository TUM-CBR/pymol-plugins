from .core.Context import Context

def __init_plugin__(app=None):

    from PyQt5.QtWidgets import QMainWindow

    from .main.Main import Main
    from .schema.SchemaContext import SchemaContext
    from .schema.visual.SchemaInstanceWidget import SchemaInstanceWidget

    context = Context()

    mw = context.run_widget(Main)

    from pymol.plugins import addmenuitemqt
    addmenuitemqt('CBR Tools', mw.show)
