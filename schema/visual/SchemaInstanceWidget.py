from PyQt5.QtWidgets import QWidget
import tempfile

from ..SchemaTaskManager import SchemaTaskManager
from .SchemaRunnerWidget import SchemaRunnerWidget
from .SchemaSelectorWidget import SchemaSelectorWidget
from .Ui_SchemaInstanceWidget import Ui_SchemaInstanceWidget

class SchemaInstanceWidget(QWidget):

    def __init__(self, schema_context, *args, **kwargs):
        super(SchemaInstanceWidget, self).__init__(*args, **kwargs)
        self.__schema_context = schema_context
        self.__ui = Ui_SchemaInstanceWidget()
        self.__ui.setupUi(self)
        self.__manager = SchemaTaskManager(schema_context, 'dummy', tempfile.gettempdir())

        self.__schema_selector = SchemaSelectorWidget(schema_context, self.__manager)
        self.__schema_runner = SchemaRunnerWidget(schema_context, self.__manager)
        self.__ui.schemaTools.addTab(self.__schema_selector, "Selector")
        self.__ui.schemaTools.addTab(self.__schema_runner, "Runner")
