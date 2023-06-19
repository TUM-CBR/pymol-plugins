from ..core.Context import Context
from .visual.SchemaInstanceWidget import SchemaInstanceWidget

def schema_raspp(context: Context):
    context.run_widget(SchemaInstanceWidget).show()