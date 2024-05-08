from ..core.Context import Context
from .visual.SequenceSearchWidget import SequenceSearchWidget

def sequence_search_main(ctx : Context):
    ctx.run_widget(SequenceSearchWidget).show()