from ..core.Context import Context
from .visual.FrankenProtBuilder import FrankenProtBuilder

def frankenprot(ctx: Context):
    ctx.run_widget(FrankenProtBuilder).show()