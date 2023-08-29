from ..core.Context import Context
from .visual.PrimerDesign import PrimerDesign

def msa_viewer(ctx : Context):
    ctx.run_widget(PrimerDesign).show()