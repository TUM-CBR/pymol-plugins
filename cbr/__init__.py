from .core.Context import Context

def __init_plugin__(app=None):

    from pymol.plugins import addmenuitemqt
    from .main.Main import Main

    context = Context()
    mw = context.run_widget(Main)
    addmenuitemqt('CBR Tools', mw.show)
