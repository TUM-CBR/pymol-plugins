from .core.Context import Context

def add_init_action(action):
    global init_actions
    init_actions.append(action)

def __init_plugin__(app=None):
    global init_actions
    from pymol.plugins import addmenuitemqt
    from .main.Main import Main

    context = Context()
    mw = context.run_widget(Main)
    addmenuitemqt('CBR Tools', mw.show)

    from .support.testing import UiRunner
    print("init cbr module")
    UiRunner.init_singleton()
