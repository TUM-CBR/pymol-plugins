from os import path
import pymol
from PyQt5.QtCore import pyqtSlot, QRegExp
from PyQt5.QtGui import QRegExpValidator
from PyQt5.QtWidgets import QWidget
from tempfile import TemporaryDirectory
from typing import List

from cbr.core.TaskManager import TaskManager

from ...clustal import Clustal
from ...core.Context import Context
from ...core.pymol import structure
from ...core import visual
from ...support.visual import as_fasta_selector
from ..raspp import schemacontacts
from ..raspp import schemaenergy
from .SchemaEnergyViewer import SchemaEnergyViewer
from .Ui_SchemaEnergyRunner import Ui_SchemaEnergyRunner

class SchemaEnergyRunner(QWidget):

    XoValidator = QRegExpValidator(QRegExp("^(\\s*\\d+\\s*(,\\s*\\d+)*)?$"))

    def __init__(self, context : Context, *args, **kwargs):
        super(SchemaEnergyRunner, self).__init__(*args, **kwargs)

        self.__ui = Ui_SchemaEnergyRunner()
        self.__ui.setupUi(self)
        visual.as_structure_selector(
            self.__ui.structuresCombo,
            self.__ui.refreshStructuresButton)

        self.__fasta_selector = as_fasta_selector(
            self.__ui.fastaTextEdit,
            self.__ui.structureSequenceCombo)

        self.__ui.shufflingPointsEdit.setValidator(SchemaEnergyRunner.XoValidator)

        self.__working_directory = TemporaryDirectory()
        self.__clustal = Clustal.get_clustal_from_context(context)
        self.__task_manager = TaskManager.from_context(context)
        self.__context = context
        self.__stop_progress_bar()
        self.__ui.runSchemaEnergyButton.clicked.connect(self.on_runSchemaEnergyButton_clicked)

    def __stop_progress_bar(self):
        self.__ui.schemaProgress.setVisible(False)
        self.__ui.runSchemaEnergyButton.setEnabled(True)

    def __start_progress_bar(self):
        self.__ui.runSchemaEnergyButton.setEnabled(False)
        self.__ui.schemaProgress.reset()
        self.__ui.schemaProgress.setMinimum(0)
        self.__ui.schemaProgress.setMaximum(0)
        self.__ui.schemaProgress.setVisible(True)

    def __del__(self):
        self.__working_directory.cleanup()

    def __structure_file(self, structure_name : str, chain_name : str) -> str:
        return path.join(
            self.__working_directory.name,
            "%s_%s.pdb" % (structure_name, chain_name)
        )

    @property
    def __msa_file(self) -> str:
        return path.join(
            self.__working_directory.name,
            "parents.aln"
        )

    @property
    def __structure_msa_file(self) -> str:
        return path.join(
            self.__working_directory.name,
            "structure.aln"
        )

    def __save_pdb(self, structure_name : str, chain_name : str) -> str:
        file_name = self.__structure_file(structure_name, chain_name)
        pymol.cmd.save(
            file_name,
            "(model %s) & (chain %s)" % (structure_name, chain_name)
        )
        return file_name

    @property
    def __contacts_file(self) -> str:
        return path.join(
            self.__working_directory.name,
            "schema_contacts.txt"
        )

    @property
    def __corssovers_file(self) -> str:
        return path.join(
            self.__working_directory.name,
            "schema_crossovers.txt"
        )

    @property
    def __schema_energy_file(self) -> str:
        return path.join(
            self.__working_directory.name,
            "schema_energy.txt"
        )

    def __save_crossovers(self, crossovers : List[int]):
        with open(self.__corssovers_file, 'w') as xo_file:
            xo_file.write(" ".join(map(str, crossovers)))
        return self.__corssovers_file

    @pyqtSlot()
    def on_runSchemaEnergyButton_clicked(self):
        (structure_name, chain_name) = self.__ui.structuresCombo.currentData()
        crossovers = [ \
            int(xo) \
            for xo in self.__ui.shufflingPointsEdit.text().split(",") \
        ]

        def task():
            return self.__run_schema_energy(
                structure_name,
                chain_name,
                crossovers)

        result = self.__task_manager.run_task(
            "schema-energy/%s/%s" % (structure_name, chain_name),
            task)

        result.on_started(self.__start_progress_bar)
        result.on_completed(self.__stop_progress_bar)
        result.on_completed(
            lambda: self.__show_results(structure_name, chain_name)
        )

    def __show_results(self, structure_name, chain_name):

        self.__context.run_widget(
            lambda _: \
                SchemaEnergyViewer(
                    self.__context,
                    structure_name,
                    chain_name,
                    self.__schema_energy_file,
                    self.__contacts_file,
                )
        ).show()

    def __run_schema_energy(
        self,
        structure_name : str,
        chain_name : str,
        crossovers : List[int]):
        
        pdb_file = self.__save_pdb(structure_name, chain_name)
        sequences = dict(self.__fasta_selector.get_items())
        self.__clustal.run_msa(
            sequences.items(),
            self.__msa_file
        )
        self.__clustal.run_msa(
            [ self.__fasta_selector.selected_sequence()
            , (structure_name, structure.get_pdb_sequence(structure_name, chain_name))
            ],
            self.__structure_msa_file
        )
        schemacontacts.main_impl({
            schemacontacts.ARG_PDB_FILE: pdb_file,
            schemacontacts.ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE: self.__msa_file,
            schemacontacts.ARG_PDB_ALIGNMENT_FILE: self.__structure_msa_file,
            schemacontacts.ARG_OUTPUT_FILE: self.__contacts_file
        })

        xo_file = self.__save_crossovers(crossovers)
        schemaenergy.main_impl({
            schemaenergy.ARG_CONTACT_FILE: self.__contacts_file,
            schemaenergy.ARG_CROSSOVER_FILE: xo_file,
            schemaenergy.ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE: self.__msa_file,
            schemaenergy.ARG_OUTPUT_FILE: self.__schema_energy_file,
            schemaenergy.ARG_PDB_ALIGNMENT_FILE: self.__structure_msa_file,
            schemaenergy.ARG_PRINT_M: True,
            schemaenergy.ARG_PRINT_E: True
        })

