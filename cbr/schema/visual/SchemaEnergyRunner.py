from Bio import AlignIO
from Bio.PDB import PDBParser
from os import path
import pymol
from PyQt5.QtCore import pyqtSlot, QRegExp
from PyQt5.QtGui import QRegExpValidator
from PyQt5.QtWidgets import QMessageBox, QWidget
from tempfile import TemporaryDirectory
from typing import Any, Dict, List, Optional

from ...core import sequence
from ...core import visual
from ...core.TaskManager import TaskManager
from ...clustal import Clustal
from ...core.Context import Context
from ...core.Qt.QtWidgets import with_error_handler
from ...support.visual import as_fasta_selector

from ..blosum import BlosumMatrix, write_matrix
from ..raspp import schemacontacts
from ..raspp import schemaenergy

from .energy import EnergySelector
from .SchemaEnergyViewer import SchemaEnergyViewer
from .SequencesPositionEditor import FragmentSelection, SequencesPositionEditor
from .substitution import BLOSUM_MATRIXES
from .Ui_SchemaEnergyRunner import Ui_SchemaEnergyRunner

class SchemaEnergyRunner(QWidget):

    XoValidator = QRegExpValidator(QRegExp("^(\\s*\\d+\\s*(,\\s*\\d+)*)?$"))

    def __init__(self, context : Context, *args: Any, **kwargs: Any):
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

        self.__fasta_selector.sequences_changed.connect(self.__on_sequences_changed)

        self.__clustal = Clustal.get_clustal_from_context(context)
        self.__task_manager = TaskManager.from_context(context)
        self.__context = context
        self.__stop_progress_bar()
        self.__ui.runSchemaEnergyButton.clicked.connect(self.on_runSchemaEnergyButton_clicked)
        self.__energy_selector = EnergySelector(self.__ui.energyScoringCombo)
        self.__positions_editor = SequencesPositionEditor(
            self.__ui.shufflingPointsTable,
            self.__ui.sequencesTable,
            self.__ui.shufflingPointBox
        )

    def __on_sequences_changed(self):
        self.__positions_editor.set_sequences(self.__fasta_selector.get_items_meta())

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

    def __msa_file(self, base_path: str) -> str:
        return path.join(
            base_path,
            "parents.aln"
        )

    def __structure_msa_file(self, base_path: str) -> str:
        return path.join(
            base_path,
            "structure.aln"
        )

    def __save_pdb(self, base_path : str, selection : visual.StructureSelection) -> str:
        file_name = self.__structure_file(base_path, selection.structure_name, selection.chain_name or "")
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

    @pyqtSlot(name="on_runSchemaEnergyButton_clicked")
    @with_error_handler()
    def on_runSchemaEnergyButton_clicked(self):
        selection = self.__structure_selector.currentSelection

        if not selection:
            raise ValueError("Select a structure!")

        fragments = self.__positions_editor.fragments()
        if fragments is None:
            raise ValueError("No fragments have been selected")

        results_directory = TemporaryDirectory()

        def task():
            return self.__run_schema_energy(
                selection,
                fragments,
                results_directory.name,
                BLOSUM_MATRIXES
            )

        result = self.__task_manager.run_task(
            "schema-energy/%s/%s" % (selection.structure_name, selection.chain_name),
            task)

        result.on_started(self.__start_progress_bar)
        result.on_completed(self.__stop_progress_bar)

        def __on_task_completed():
            if result.error:
                QMessageBox.critical(self, "Error", str(result.error))
            else:
                self.__show_results(selection, results_directory)
        result.on_completed(__on_task_completed)

    def __show_results(
        self,
        structure_seleciton : visual.StructureSelection,
        results_folder : TemporaryDirectory[Any]
    ):

        self.__context.run_widget(
            lambda _: \
                SchemaEnergyViewer(
                    self.__context,
                    structure_seleciton,
                    self.__schema_energy_file(results_folder.name),
                    self.__contacts_file(results_folder.name),
                    self.__msa_file(results_folder.name),
                    results_folder
                )
        ).show()

    def __with_blosum(self, base_path : str, matrix : Dict[Any, Any]) -> str:
        location = path.join(base_path, "blosum.json")
        with open(location, 'w') as blosum:
            write_matrix(matrix, blosum)

        return location

    def __build_msa_from_fragmetns(
        self,
        fragments: FragmentSelection,
        out_file: str
    ) -> List[int]:
        alignment = self.__clustal.run_msa_fragments(
            fragments.fragments,
            result_order = fragments.order
        )

        with open(out_file, 'w') as handle:
            AlignIO.write(
                alignment.alignment,
                handle,
                format='clustal'
            )
        return alignment.fragments_positions

    def __run_schema_energy(
        self,
        structure_selection : visual.StructureSelection,
        fragments: FragmentSelection,
        base_path : str,
        blosum : Optional[Dict[str, BlosumMatrix]]
    ):        
        pdb_file = self.__save_pdb(base_path, structure_selection)

        parser = PDBParser()
        pdb_structure = parser.get_structure(
            structure_selection.structure_name,
            pdb_file
        )
        residues : Any = pdb_structure.get_residues()
        pdb_seq = "".join(sequence.residue_to_1(res.resname) for res in residues)

        parents_msa = self.__msa_file(base_path)
        structure_msa = self.__structure_msa_file(base_path)

        crossovers = self.__build_msa_from_fragmetns(
            fragments,
            parents_msa
        )

        self.__clustal.run_msa(
            [ self.__fasta_selector.selected_sequence()
            , (structure_selection.structure_name, pdb_seq)
            ],
            structure_msa
        )
        schemacontacts.main_impl({
            schemacontacts.ARG_PDB_FILE: pdb_file,
            schemacontacts.ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE: parents_msa,
            schemacontacts.ARG_PDB_ALIGNMENT_FILE: structure_msa,
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
            schemaenergy_args[schemaenergy.ARG_BLOSUM] = self.__with_blosum(base_path, blosum)

        schemaenergy.main_impl(schemaenergy_args)

