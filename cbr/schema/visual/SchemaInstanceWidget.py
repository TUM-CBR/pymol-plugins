from PyQt5.QtWidgets import QWidget

from ...core.Context import Context
from ..SchemaContext import SchemaContext
from ..SchemaTaskManager import SchemaTaskManager
from .SchemaRunnerWidget import SchemaRunnerWidget
from .SchemaSelectorWidget import SchemaSelectorWidget
from .Ui_SchemaInstanceWidget import Ui_SchemaInstanceWidget

def create_context():
    return SchemaContext()

K_CONTEXT = 'schema-context'

class SchemaInstanceWidget(QWidget):

    def __init__(self, context : Context, *args, **kwargs):
        super(SchemaInstanceWidget, self).__init__(*args, **kwargs)
        schema_context = context.create_or_load(K_CONTEXT, create_context)
        self.__ui = Ui_SchemaInstanceWidget()
        self.__ui.setupUi(self)
        self.__manager = SchemaTaskManager(schema_context)

        self.__schema_selector = SchemaSelectorWidget(schema_context, self.__manager)
        self.__schema_runner = SchemaRunnerWidget(schema_context, self.__manager)
        self.__ui.schemaTools.addTab(self.__schema_selector, "Selector")
        self.__ui.schemaTools.addTab(self.__schema_runner, "Runner")
