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

from .SchemaRunnerUI import Ui_SchemaRunner
from ..SchemaTaskManager import SchemaTaskManager

class SchemaRunner(QWidget):

    def __init__(self, schema_context, *args, **kwargs):
        super(SchemaRunner, self).__init__(*args, **kwargs)
        self.__schema_context = schema_context
        self.__ui = Ui_SchemaRunner()
        self.__ui.setupUi(self)
        self.__manager = SchemaTaskManager(schema_context, 'dummy', tempfile.gettempdir())

    def __refresh_structures(self):
        self.__ui.structuresCombo.clear()
        self.__ui.structuresCombo.addItems(pymol.cmd.get_names())

    def on_refreshButton_clicked(self):
        self.__refresh_structures()

    def on_runSchemaButton_clicked(self):

        try:
            pymol.cmd.save(self.__manager.pdb_file, self.__ui.structuresCombo.currentText())
            pymol.cmd.load(self.__manager.pdb_file)

            sequences = "%s\n\n>%s\n%s" % (self.__ui.sequencesText.toPlainText(), self.__manager.structure_name, pymol.cmd.get_fastastr(self.__manager.structure_name))
            
            self.__manager.run_schema(sequences, 5, 30)
        except QuietException:
            self.__schema_context.raise_error_message("The structure '%s' is invalid" % self.__ui.structuresCombo.currentText)
        
