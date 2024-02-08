import json
from pathlib import Path
from PyQt5.QtCore import pyqtSlot, QUrl
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtWidgets import QWidget
from os import path
from pymol import cmd
from typing import Iterable, Optional

from ...core.executable import Executable, ExecutableGroupResult, ExecutableProcess, ExecutableProcessGroup
from ...core.Context import Context
from ...core.Qt.QtWidgets import show_error
from ..data import MpnnSpec, mpnn_selection
from ..dispatch import parse_multiple_chains
from .Ui_MpnnViewer import Ui_MpnnViewer

class MpnnViewer(QWidget):

    def __init__(
        self,
        context: Context,
        mpnn: Executable,
        spec: MpnnSpec
    ) -> None:
        super().__init__()
        self.__context = context
        self.__working_directory = self.__context.create_temporary_directory()
        self.__mpnn = mpnn
        self.__spec = spec
        self.__ui = Ui_MpnnViewer()
        self.__ui.setupUi(self)
        self.__processess : Optional[ExecutableProcessGroup] = None
        self.__ui.resultsFolderButton.clicked.connect(self.__on_results_folder_clicked)

        self.__init_widget()

    @pyqtSlot()
    def __on_results_folder_clicked(self):
        QDesktopServices.openUrl(
            QUrl.fromLocalFile(self.__working_directory)
        )

    def __get_chains_jonsl_path(self) -> str:
        return path.join(self.__working_directory, "chains.jsonl")
    
    def __get_assigned_chains_jonsl_path(self) -> str:
        return path.join(self.__working_directory, "assigned.jsonl")
    
    def __get_fixed_positions_jonsl_path(self) -> str:
        return path.join(self.__working_directory, "fixed.jsonl")
    
    def __get_model_location(self, model: str) -> str:
        return path.join(self.__working_directory, f"{model}.pdb")
    
    def __get_results_location(self) -> str:
        return path.join(self.__working_directory, "results")

    def __save_models(self, models: Iterable[str]):

        for model in models:
            cmd.save(
                self.__get_model_location(model),
                mpnn_selection(model)
            )

        parse_multiple_chains.main(
            folder_with_pdbs_path=self.__working_directory,
            save_path=self.__get_chains_jonsl_path(),
            ca_only=False
        )

    def __save_fixed_chains(self):
        spec = self.__spec
        with open(self.__get_assigned_chains_jonsl_path(), 'w') as jonsl:
            json.dump(spec.get_chains_jsonl(), jonsl)

    def __save_fixed_positions(self):
        spec = self.__spec
        with open(self.__get_fixed_positions_jonsl_path(), 'w') as jsonl:
            json.dump(spec.get_positions_jsonl(), jsonl)

    def __populate_fasta(self):

        result = ""
        for fasta in Path(self.__get_results_location()).rglob("*.[Ff][aA]"):
            with fasta.open('r') as f:
                result += f.read() + "\n\n"

        self.__ui.sequencesTextEdit.setPlainText(result)

    @pyqtSlot(object)
    def __on_complete(self, result: ExecutableGroupResult):

        processes = self.__processess

        if processes is None:
            return

        processes.on_complete.disconnect(self.__on_complete)
        processes.close()
        processes = None

        self.__ui.progressBar.hide()

        errors = result.get_errors()
        if errors is not None:
            show_error(self, "Protein MPNN Failed", errors)

        self.__populate_fasta()

    def __run_mpnn(self):
        
        self.__processess = ExecutableProcessGroup.create(
            ExecutableProcess.create_process(
                self.__mpnn,
                [
                    "--pdb_path", self.__get_model_location(model),
                    "--out_folder", self.__get_results_location(),
                    "--jsonl_path", self.__get_chains_jonsl_path(),
                    "--chain_id_jsonl", self.__get_assigned_chains_jonsl_path(),
                    "--fixed_positions_jsonl", self.__get_fixed_positions_jonsl_path(),
                    "--num_seq_per_target", str(self.__spec.num_seqs)
                ]
            )
            for model in self.__spec.get_models()
        )
        self.__processess.on_complete.connect(self.__on_complete)
        self.__processess.start()

    def __init_widget(self):
        self.__save_models(self.__spec.get_models())
        self.__save_fixed_chains()
        self.__save_fixed_positions()
        self.__run_mpnn()




