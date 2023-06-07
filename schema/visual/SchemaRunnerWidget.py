# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'SchemaRunner.ui',
# licensing of 'SchemaRunner.ui' applies.
#
# Created: Mon Jun  5 09:16:29 2023
#      by: pyside2-uic  running on PySide2 5.13.2
#
# WARNING! All changes made in this file will be lost!

import pymol
from pymol.parsing import QuietException
from PyQt5.QtWidgets import QWidget
import tempfile

from .Ui_SchemaRunnerWidget import Ui_SchemaRunnerWidget
from ..SchemaTaskManager import SchemaTaskManager

input_validation_error = \
    """The values of the input fields could not be read:
    * Please ensure that the crossovers are a comma separated integers
    """

class SchemaRunnerWidget(QWidget):

    SEQUENCES_KEY = 'sequences'

    def __init__(self, schema_context, manager, *args, **kwargs):
        super(SchemaRunnerWidget, self).__init__(*args, **kwargs)
        self.__schema_context = schema_context
        self.__ui = Ui_SchemaRunnerWidget()
        self.__ui.setupUi(self)
        self.__manager = manager
        self.__manager.is_busy_signal.connect(self.__on_busy_status_changed)
        self.__load_previous_sequences()
        self.__on_busy_status_changed(False)

    def __on_busy_status_changed(self, is_busy : bool):
        if is_busy:
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

    def __refresh_structures(self):
        self.__ui.structuresCombo.clear()
        self.__ui.structuresCombo.addItems(pymol.cmd.get_names())

    def on_refreshButton_clicked(self):
        self.__refresh_structures()

    def __get_pdb_sequence(self):
        result = []
        pymol.cmd.iterate(
            "(%s) & guide & alt +A" % self.__manager.structure_name,
            'result.append(resn)',
            space={'result': result}
        )
        legend = pymol.exporting._resn_to_aa

        return "".join(legend[resn] for resn in result)

    def __validate_crossovers(self):
        xos = self.__ui.crossoversText.text()
        return [int(x) for x in xos.split(",")]


    def on_runSchemaButton_clicked(self):

        self.__save_sequences()

        try:
            xos = self.__validate_crossovers()
            pymol.cmd.save(self.__manager.pdb_file, self.__ui.structuresCombo.currentText())
            pymol.cmd.load(self.__manager.pdb_file)

            sequences = "%s\n\n>%s\n%s" % (self.__ui.sequencesText.toPlainText(), self.__manager.structure_name, self.__get_pdb_sequence())
            
            self.__manager.run_schema(sequences, xos, 3)
        except QuietException:
            self.__schema_context.raise_error_message("The structure '%s' is invalid" % self.__ui.structuresCombo.currentText)
        except ValueError:
            self.__schema_context.raise_error_message(input_validation_error)
        
