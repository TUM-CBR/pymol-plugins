from ..core.Context import Context
from .visual.CoevolutionRunner import CoevolutionRunner

def coevolution_runner(ctx : Context):
    ctx.run_widget(CoevolutionRunner).show()