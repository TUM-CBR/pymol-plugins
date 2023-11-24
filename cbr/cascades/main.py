from ..core.Context import Context
from .visual.CascadesMain import CascadesMain

def cascades_main(ctx : Context):
    ctx.run_widget(CascadesMain).show()