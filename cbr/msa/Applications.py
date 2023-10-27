from ..core.Context import Context
from .visual.MsaViewer import MsaViewer
from .visual.MsaCleaner import MsaCleaner

def msa_viewer(ctx : Context):
    ctx.run_widget(MsaViewer).show()

def msa_cleaner(ctx : Context):
    ctx.run_widget(MsaCleaner).show()