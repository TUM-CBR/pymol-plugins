from ..core.Context import Context
from .visual.Dssp import Dssp

def dssp_runner(ctx : Context):
    ctx.run_widget(
        lambda _: Dssp()
    ).show()