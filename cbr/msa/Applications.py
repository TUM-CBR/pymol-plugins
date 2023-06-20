from ..core.Context import Context
from .visual.MsaViewer import MsaViewer

def msa_viewer(ctx : Context):
    ctx.run_widget(MsaViewer).show()