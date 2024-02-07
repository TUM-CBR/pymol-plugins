import json
from PyQt5.QtWidgets import QWidget
from os import path
from pymol import cmd
from typing import Iterable

from ...core.executable import Executable, ExecutableProcess
from ...core.Context import Context
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
        self.__running = []

        self.__init_widget()

    def __get_chains_jonsl_path(self) -> str:
        return path.join(self.__working_directory, "chains.jsonl")
    
    def __get_assigned_chains_jonsl_path(self) -> str:
        return path.join(self.__working_directory, "assigned.jsonl")
    
    def __get_fixed_positions_jonsl_path(self) -> str:
        return path.join(self.__working_directory, "fixed.jsonl")
    
    def __get_model_location(self, model: str) -> str:
        return path.join(self.__working_directory, f"{model}.pdb")
    
    def __get_results_location(self) -> str:
        return path.join(self.__working_directory, "out")

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

    def __run_mpnn(self):
        self.__running = processess = [
            ExecutableProcess.create_process(
                self.__mpnn,
                [
                    "--pdb_path", self.__get_model_location(model),
                    "--out_folder", self.__get_results_location(),
                    "--jsonl_path", self.__get_chains_jonsl_path(),
                    "--chain_id_jsonl", self.__get_assigned_chains_jonsl_path(),
                    "--fixed_positions_jsonl", self.__get_fixed_positions_jonsl_path()
                ]
            )
            for model in self.__spec.get_models()
        ]

    def __init_widget(self):
        self.__save_models(self.__spec.get_models())
        self.__save_fixed_chains()
        self.__save_fixed_positions()




