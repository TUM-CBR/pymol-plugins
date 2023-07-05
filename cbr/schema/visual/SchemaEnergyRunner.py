from os import path
import pymol
from PyQt5.QtCore import pyqtSlot, QRegExp
from PyQt5.QtGui import QRegExpValidator
from PyQt5.QtWidgets import QWidget
from tempfile import TemporaryDirectory
from typing import List, Optional

from ...core.TaskManager import TaskManager
from ...clustal import Clustal
from ...core.Context import Context
from ...core.pymol import structure
from ...core import visual
from ...support.visual import as_fasta_selector
from ..raspp import contacts
from ..raspp import schemacontacts
from ..raspp import schemaenergy
from .SchemaEnergyViewer import SchemaEnergyViewer
from .Ui_SchemaEnergyRunner import Ui_SchemaEnergyRunner

interactions = {
    'SCHEMA classic' : None,
    'Simplified Pysics' : [
        contacts.electrostatic_interactions,
        contacts.van_der_waals
    ]
}

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

        self.__clustal = Clustal.get_clustal_from_context(context)
        self.__task_manager = TaskManager.from_context(context)
        self.__context = context
        self.__stop_progress_bar()
        self.__ui.runSchemaEnergyButton.clicked.connect(self.on_runSchemaEnergyButton_clicked)

        for (i, (name, data)) in enumerate(interactions.items()):
            self.__ui.energyScoringCombo.insertItem(
                i,
                name,
                userData=data
            )

    def __stop_progress_bar(self):
        self.__ui.schemaProgress.setVisible(False)
        self.__ui.runSchemaEnergyButton.setEnabled(True)

    def __start_progress_bar(self):
        self.__ui.runSchemaEnergyButton.setEnabled(False)
        self.__ui.schemaProgress.reset()
        self.__ui.schemaProgress.setMinimum(0)
        self.__ui.schemaProgress.setMaximum(0)
        self.__ui.schemaProgress.setVisible(True)

    def __structure_file(self, base_path : str, structure_name : str, chain_name : str) -> str:
        return path.join(
            base_path,
            "%s_%s.pdb" % (structure_name, chain_name)
        )

    def __msa_file(self, base_path) -> str:
        return path.join(
            base_path,
            "parents.aln"
        )

    def __structure_msa_file(self, base_path) -> str:
        return path.join(
            base_path,
            "structure.aln"
        )

    def __save_pdb(self, base_path : str, structure_name : str, chain_name : str) -> str:
        file_name = self.__structure_file(base_path, structure_name, chain_name)
        pymol.cmd.save(
            file_name,
            "(model %s) & (chain %s)" % (structure_name, chain_name)
        )
        return file_name

    def __contacts_file(self, base_path : str) -> str:
        return path.join(
            base_path,
            "schema_contacts.txt"
        )

    def __corssovers_file(self, base_path : str) -> str:
        return path.join(
            base_path,
            "schema_crossovers.txt"
        )

    def __schema_energy_file(self, base_path : str) -> str:
        return path.join(
            base_path,
            "schema_energy.txt"
        )

    def __schema_interactions_file(self, base_path : str) -> str:
        return path.join(
            base_path,
            "schema_interactions.json"
        )

    def __save_crossovers(self, base_path : str, crossovers : List[int]) -> str:
        with open(self.__corssovers_file(base_path), 'w') as xo_file:
            xo_file.write(" ".join(map(str, crossovers)))
        return self.__corssovers_file(base_path)

    @pyqtSlot()
    def on_runSchemaEnergyButton_clicked(self):
        (structure_name, chain_name) = self.__ui.structuresCombo.currentData()
        crossovers = [ \
            int(xo) \
            for xo in self.__ui.shufflingPointsEdit.text().split(",") \
        ]

        results_directory = TemporaryDirectory()

        def task():
            return self.__run_schema_energy(
                structure_name,
                chain_name,
                crossovers,
                results_directory.name
            )

        result = self.__task_manager.run_task(
            "schema-energy/%s/%s" % (structure_name, chain_name),
            task)

        result.on_started(self.__start_progress_bar)
        result.on_completed(self.__stop_progress_bar)
        result.on_completed(
            lambda: self.__show_results(structure_name, chain_name, results_directory)
        )

    def __show_results(
        self,
        structure_name : str,
        chain_name : str,
        results_folder : TemporaryDirectory
        ):

        self.__context.run_widget(
            lambda _: \
                SchemaEnergyViewer(
                    self.__context,
                    structure_name,
                    chain_name,
                    self.__schema_energy_file(results_folder.name),
                    self.__contacts_file(results_folder.name),
                    results_folder
                )
        ).show()

    def __get_interactions(self, base_path : str) -> Optional[str]:

        scoring = self.__ui.energyScoringCombo.currentData()

        if scoring is None:
            return None
        else:
            f_name = self.__schema_interactions_file(base_path)
            with open(f_name, 'w') as f_interactions:
                contacts.write_interactions(scoring, f_interactions)
            return f_name

    def __run_schema_energy(
        self,
        structure_name : str,
        chain_name : str,
        crossovers : List[int],
        base_path : str):
        
        pdb_file = self.__save_pdb(base_path, structure_name, chain_name)
        sequences = dict(self.__fasta_selector.get_items())
        self.__clustal.run_msa(
            sequences.items(),
            self.__msa_file(base_path)
        )
        self.__clustal.run_msa(
            [ self.__fasta_selector.selected_sequence()
            , (structure_name, structure.get_pdb_sequence(structure_name, chain_name))
            ],
            self.__structure_msa_file(base_path)
        )
        schemacontacts.main_impl({
            schemacontacts.ARG_PDB_FILE: pdb_file,
            schemacontacts.ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE: self.__msa_file(base_path),
            schemacontacts.ARG_PDB_ALIGNMENT_FILE: self.__structure_msa_file(base_path),
            schemacontacts.ARG_OUTPUT_FILE: self.__contacts_file(base_path),
            schemacontacts.ARG_INTERACTIONS: self.__get_interactions(base_path)
        })

        xo_file = self.__save_crossovers(base_path, crossovers)
        schemaenergy.main_impl({
            schemaenergy.ARG_CONTACT_FILE: self.__contacts_file(base_path),
            schemaenergy.ARG_CROSSOVER_FILE: xo_file,
            schemaenergy.ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE: self.__msa_file(base_path),
            schemaenergy.ARG_OUTPUT_FILE: self.__schema_energy_file(base_path),
            schemaenergy.ARG_PDB_ALIGNMENT_FILE: self.__structure_msa_file(base_path),
            schemaenergy.ARG_PRINT_M: True,
            schemaenergy.ARG_PRINT_E: True
        })

