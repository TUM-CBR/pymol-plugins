from ..core.Context import Context
from .visual.ChimerasGenerator import ChimerasGenerator

def chimeras_generator(ctx : Context):
    ctx.run_widget(ChimerasGenerator).show()