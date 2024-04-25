from ..core.Context import Context
from .visual.Coevolution import Coevolution

def coevolution_runner(ctx : Context):
    ctx.run_widget(Coevolution).show()