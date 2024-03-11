from ..core.Context import Context
from .visual.CavityFinder import CavityFinder

def main(context: Context):
    context.run_widget(CavityFinder).show()