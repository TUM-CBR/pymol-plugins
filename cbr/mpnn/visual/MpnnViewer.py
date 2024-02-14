from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import json
from pathlib import Path
from PyQt5.QtCore import pyqtSlot, QUrl
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtWidgets import QWidget
from os import path
from pymol import cmd
import re
from typing import cast, Dict, Iterable, List, Optional

from ...core.executable import Executable, ExecutableGroupResult, ExecutableProcess, ExecutableProcessGroup
from ...core.Context import Context
from ...core.Qt.QtWidgets import show_error
from ...core.pymol.structure import StructureSelection
from ..data import MpnnSpec, mpnn_selection
from ..dispatch import parse_multiple_chains
from .Ui_MpnnViewer import Ui_MpnnViewer

DESIGNED_CHAINS_RE = re.compile(r"designed_chains=\[[a-zA-z',]+\]")

def clean_id(id: Optional[str]) -> Optional[str]:

    if id is None:
        return None
    
    return id.replace(",", "").upper()

def get_selection(record: SeqRecord) -> Optional[StructureSelection]:
    model = clean_id(record.id)
    chains = DESIGNED_CHAINS_RE.findall(record.description)
    models: List[str] = cmd.get_names()
    model_ix = next(
        (ix for ix, candidate in enumerate(models) if candidate.upper() == model),
        None
    )

    if model is None \
        or len(chains) < 1 \
        or model_ix is None:
        return None
    
    model = models[model_ix]
    seq_chains = [chain.replace("'", "").upper() for chain in chains[0]]
    available_chains = [chain.upper() for chain in cmd.get_chains(model)]

    for chain in seq_chains:

        if chain in available_chains:
            return StructureSelection(
                structure_name=model,
                chain_name=chain,
                segment_identifier=None
            )
        
    return None

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
    
    def __get_excluded_residues_jonsl_path(self) -> str:
        return path.join(self.__working_directory, "excluded.jsonl")
    
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

        seqs: List[SeqRecord] = []
        selections: Dict[int, StructureSelection] = {}
        for fasta in Path(self.__get_results_location()).rglob("*.[Ff][aA]"):
            current = list(cast(Iterable[SeqRecord], SeqIO.parse(fasta, format='fasta')))

            if len(current) < 0:
                continue

            main = current[0]
            start_ix = len(seqs)
            seqs += current
            
            # Find out which of the currently loaded structures was used
            # as a template to generate the given sequences
            selection = get_selection(main)
            if selection is None:
                continue

            for i in range(start_ix, len(seqs)):
                selections[i] = selection


        self.__ui.fastaViewer.set_sequences(seqs, structure_mappings=selections)

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

    def __get_excluded_args(self) -> List[str]:
        excluded = self.__spec.get_excluded_jonsl()
        if len(excluded) < 1:
            return []
        
        excluded_file = self.__get_excluded_residues_jonsl_path()
        with open(excluded_file, 'w') as out_stream:
            json.dump(excluded, out_stream)

        return ["--omit_AA_jsonl", excluded_file]

    def __run_mpnn(self):

        args = self.__spec.mpnn_args
        self.__processess = ExecutableProcessGroup.create(
            ExecutableProcess.create_process(
                self.__mpnn,
                [
                    "--pdb_path", self.__get_model_location(model),
                    "--out_folder", self.__get_results_location(),
                    "--jsonl_path", self.__get_chains_jonsl_path(),
                    "--chain_id_jsonl", self.__get_assigned_chains_jonsl_path(),
                    "--fixed_positions_jsonl", self.__get_fixed_positions_jonsl_path(),
                    "--num_seq_per_target", str(self.__spec.num_seqs),
                    "--backbone_noise", str(args.backbone_noise),
                    "--sampling_temp", str(args.sampling_temperature)
                ] \
                + self.__get_excluded_args() \
                + (["--use_soluble_model"] if args.use_soluble_model else [])
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




