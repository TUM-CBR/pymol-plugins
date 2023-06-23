from ..core.Context import Context
from .visual.SchemaEnergyRunner import SchemaEnergyRunner
from .visual.SchemaInstanceWidget import SchemaInstanceWidget

def schema_raspp(context: Context):
    context.run_widget(SchemaInstanceWidget).show()

def schema_energy(context: Context):
    context.run_widget(SchemaEnergyRunner).show()