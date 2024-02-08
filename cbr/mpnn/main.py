from ..core.Context import Context
from .visual.MpnnRunner import MpnnRunner

def mpnn_main(ctx : Context):
    ctx.run_widget(MpnnRunner).show()