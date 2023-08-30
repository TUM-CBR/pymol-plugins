from ..core.Context import Context
from .visual.PrimerDesign import PrimerDesign

def primer_generator(ctx : Context):
    ctx.run_widget(
        lambda _: PrimerDesign()
    ).show()