from ..core.Context import Context
from .visual.FastaViewerApp import FastaViewerApp
from .visual.MsaViewer import MsaViewer
from .visual.MsaCleaner import MsaCleaner

def msa_viewer(ctx : Context):
    ctx.run_widget(MsaViewer).show()

def msa_cleaner(ctx : Context):
    ctx.run_widget(MsaCleaner).show()

def fasta_viewer(ctx: Context):
    ctx.run_widget(lambda _ctx: FastaViewerApp()).show()