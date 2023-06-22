import datetime
import pymol
from pymol.parsing import QuietException
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QIntValidator
from PyQt5.QtWidgets import QWidget

from ...core import visual
from ...core.pymol import structure
from ..SchemaTaskManager import SchemaTask, SchemaTaskManager

from .Ui_SchemaRunnerWidget import Ui_SchemaRunnerWidget

input_validation_error = \
    """The values of the input fields could not be read:
    * Please ensure that the crossovers are a comma separated integers
    """

class SchemaRunnerWidget(QWidget):

    SEQUENCES_KEY = 'sequences'

    def __init__(self, schema_context, manager : SchemaTaskManager, *args, **kwargs):
        super(SchemaRunnerWidget, self).__init__(*args, **kwargs)
        self.__schema_context = schema_context
        self.__ui = Ui_SchemaRunnerWidget()
        self.__ui.setupUi(self)

        self.__ui.minLengthText.setValidator(QIntValidator(1,9999))
        self.__ui.maxLengthText.setValidator(QIntValidator(1,9999))

        self.__manager = manager
        self.__manager.is_busy_signal.connect(self.__on_busy_status_changed)
        self.__load_previous_sequences()
        self.__on_busy_status_changed(False)
        visual.as_structure_selector(self.__ui.structuresCombo, self.__ui.refreshButton)

    def __on_busy_status_changed(self, is_busy : int):
        if is_busy > 0:
            self.__ui.runSchemaButton.setEnabled(False)
            self.__ui.schemaProgress.setRange(0,0)
            self.__ui.schemaProgress.setVisible(True)
        else:
            self.__ui.runSchemaButton.setEnabled(True)
            self.__ui.schemaProgress.setVisible(False)

    def __load_previous_sequences(self):
        previous = self.__manager.load_resource(SchemaRunnerWidget.SEQUENCES_KEY)

        if previous:
            self.__ui.sequencesText.setPlainText(previous)

    def __save_sequences(self):
        self.__manager.save_resource(
            SchemaRunnerWidget.SEQUENCES_KEY,
            self.__ui.sequencesText.toPlainText()
        )

    def __get_pdb_sequence(self, structure_name):
        return structure.get_pdb_sequence(structure_name)

    def __validate_crossovers(self):
        xos = self.__ui.crossoversText.text()
        return [int(x) for x in xos.split(",")]

    def __get_pdb_file_name(self, name : str) -> str:
        return self.__manager.get_pdb_file_name(name)

    def __get_structure_name(self, name : str) -> str:
        return SchemaTask.get_structure_name(name)

    @pyqtSlot()
    def on_runSchemaButton_clicked(self):

        self.__save_sequences()

        try:
            name = datetime.datetime.now().strftime('%Y_%m_%d_%H%M%S')
            pdb_file = self.__get_pdb_file_name(name)
            structure_name = self.__get_structure_name(name)
            xos = self.__validate_crossovers()
            min_length = int(self.__ui.minLengthText.text())
            max_length = int(self.__ui.maxLengthText.text())
            
            pymol.cmd.save(pdb_file, "(model %s) & (chain %s)" % self.__ui.structuresCombo.currentData())
            pymol.cmd.load(pdb_file)
            sequences = "%s\n\n>%s\n%s" % (self.__ui.sequencesText.toPlainText(), structure_name, self.__get_pdb_sequence(structure_name))

            self.__manager.run_schema(name, sequences, xos, min_length, max_length)
            
        except QuietException:
            self.__schema_context.raise_error_message("The structure '%s' is invalid" % self.__ui.structuresCombo.currentText)
        except ValueError:
            self.__schema_context.raise_error_message(input_validation_error)
        
