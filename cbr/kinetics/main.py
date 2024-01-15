from .visual.KineticsOptimizer import KineticsOptimizer

from ..core.Context import Context
from ..core.Qt.QtWidgets import show_error, show_info

REQ_MESSAGE = """To use this application, we need to install 'matplotlib' (https://matplotlib.org/) in order to generate graphs.
Once you close this window, PyMOL will ask you if you like to install it. Please accept the dialog box
and wait for some time (PyMOL will freeze) until you are informed that the library has been installed.
"""

REQ_ERROR = """The requirement 'matplotlib' failed to install."""

def install_deps():
    show_info(
        None,
        "App Requirements",
        REQ_MESSAGE
    )

    from pymol.Qt import utils

    utils.conda_ask_install("matplotlib", channel="conda-forge")

    try:
        import matplotlib
    except ModuleNotFoundError:
        show_error(
            None,
            "Installation Error",
            REQ_ERROR
        )
        return False

    return True


def enzyme_kinetics(ctx: Context):

    try:
        import matplotlib
    except ModuleNotFoundError:
        if not install_deps():
            return

    from .visual.KineticsInput import KineticsInput
    ctx.run_widget(
        KineticsOptimizer
    ).show()
