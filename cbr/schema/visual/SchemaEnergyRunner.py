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

from ..blosum import BlosumMatrix, write_matrix
from ..raspp import schemacontacts
from ..raspp import schemaenergy

from .energy import EnergySelector
from .SchemaEnergyViewer import SchemaEnergyViewer
from .substitution import SubstitutionSelector
from .Ui_SchemaEnergyRunner import Ui_SchemaEnergyRunner

class SchemaEnergyRunner(QWidget):

    XoValidator = QRegExpValidator(QRegExp("^(\\s*\\d+\\s*(,\\s*\\d+)*)?$"))

    def __init__(self, context : Context, *args, **kwargs):
        super(SchemaEnergyRunner, self).__init__(*args, **kwargs)

        self.__ui = Ui_SchemaEnergyRunner()
        self.__ui.setupUi(self)
        self.__structure_selector = visual.as_structure_selector(
            self.__ui.structuresCombo,
            self.__ui.refreshStructuresButton,
            copy_button = self.__ui.copySequenceButton
        )

        self.__fasta_selector = as_fasta_selector(
            self.__ui.fastaTextEdit,
            self.__ui.structureSequenceCombo
        )

        self.__ui.shufflingPointsEdit.setValidator(SchemaEnergyRunner.XoValidator)

        self.__clustal = Clustal.get_clustal_from_context(context)
        self.__task_manager = TaskManager.from_context(context)
        self.__context = context
        self.__stop_progress_bar()
        self.__ui.runSchemaEnergyButton.clicked.connect(self.on_runSchemaEnergyButton_clicked)
        self.__energy_selector = EnergySelector(self.__ui.energyScoringCombo)
        self.__substitution_selector = SubstitutionSelector(self.__ui.substitutionSelector)

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

    def __save_pdb(self, base_path : str, selection : visual.StructureSelection) -> str:
        file_name = self.__structure_file(base_path, selection.structure_name, selection.chain_name)
        pymol.cmd.save(
            file_name,
            selection.selection
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

    def __save_crossovers(self, base_path : str, crossovers : List[int]) -> str:
        with open(self.__corssovers_file(base_path), 'w') as xo_file:
            xo_file.write(" ".join(map(str, crossovers)))
        return self.__corssovers_file(base_path)

    @pyqtSlot()
    def on_runSchemaEnergyButton_clicked(self):
        selection = self.__structure_selector.currentSelection

        if not selection:
            raise ValueError("Select a structure!")

        crossovers = [ \
            int(xo) \
            for xo in self.__ui.shufflingPointsEdit.text().split(",") \
        ]

        results_directory = TemporaryDirectory()

        def task():
            return self.__run_schema_energy(
                selection,
                crossovers,
                results_directory.name,
                self.__substitution_selector.selection
            )

        result = self.__task_manager.run_task(
            "schema-energy/%s/%s" % (selection.structure_name, selection.chain_name),
            task)

        result.on_started(self.__start_progress_bar)
        result.on_completed(self.__stop_progress_bar)
        result.on_completed(
            lambda: self.__show_results(selection, results_directory)
        )

    def __show_results(
        self,
        structure_seleciton : visual.StructureSelection,
        results_folder : TemporaryDirectory
        ):

        self.__context.run_widget(
            lambda _: \
                SchemaEnergyViewer(
                    self.__context,
                    structure_seleciton,
                    self.__schema_energy_file(results_folder.name),
                    self.__contacts_file(results_folder.name),
                    results_folder
                )
        ).show()

    def __with_blosum(self, base_path : str, matrix : dict) -> str:
        location = path.join(base_path, "blosum.json")
        with open(location, 'w') as blosum:
            write_matrix(matrix, blosum)

        return location

    def __run_schema_energy(
        self,
        structure_selection : visual.StructureSelection,
        crossovers : List[int],
        base_path : str,
        blosum : Optional[BlosumMatrix]
    ):
        
        pdb_file = self.__save_pdb(base_path, structure_selection)
        sequences = dict(self.__fasta_selector.get_items())
        self.__clustal.run_msa(
            sequences.items(),
            self.__msa_file(base_path)
        )
        self.__clustal.run_msa(
            [ self.__fasta_selector.selected_sequence()
            , (structure_selection.structure_name, structure.get_selection_sequece(structure_selection.selection))
            ],
            self.__structure_msa_file(base_path)
        )
        schemacontacts.main_impl({
            schemacontacts.ARG_PDB_FILE: pdb_file,
            schemacontacts.ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE: self.__msa_file(base_path),
            schemacontacts.ARG_PDB_ALIGNMENT_FILE: self.__structure_msa_file(base_path),
            schemacontacts.ARG_OUTPUT_FILE: self.__contacts_file(base_path),
            schemacontacts.ARG_INTERACTIONS: self.__energy_selector.write_interactions(base_path)
        })

        xo_file = self.__save_crossovers(base_path, crossovers)
        schemaenergy_args = {
            schemaenergy.ARG_CONTACT_FILE: self.__contacts_file(base_path),
            schemaenergy.ARG_CROSSOVER_FILE: xo_file,
            schemaenergy.ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE: self.__msa_file(base_path),
            schemaenergy.ARG_OUTPUT_FILE: self.__schema_energy_file(base_path),
            schemaenergy.ARG_PDB_ALIGNMENT_FILE: self.__structure_msa_file(base_path),
            schemaenergy.ARG_PRINT_M: True,
            schemaenergy.ARG_PRINT_E: True,
        }

        if blosum:
            schemaenergy_args[schemaenergy.ARG_DISRUPTION] = self.__with_blosum(base_path, blosum)

        schemaenergy.main_impl(schemaenergy_args)

