from ..core.Context import Context
from .visual.FrankenProt import FrankenProt

def frankenprot(ctx: Context):
    ctx.run_widget(FrankenProt).show()