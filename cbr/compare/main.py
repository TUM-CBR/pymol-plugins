from ..core.Context import Context
from ..support.structure.StructureCompare import StructureCompare

def main(context: Context):
    context.run_widget(StructureCompare).show()