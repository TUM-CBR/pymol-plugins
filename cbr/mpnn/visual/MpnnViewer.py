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
from typing import Any, Sequence, cast, Dict, Iterable, List, Optional

from ...core.executable import Executable, ExecutableGroupResult, ExecutableProcess, ExecutableProcessGroup
from ...core.Context import Context
from ...core.Qt.QtWidgets import show_error
from ...core.pymol.structure import StructureSelection
from ...support.fasta.visual.FastaViewer import SequenceStructures
from ..data import MpnnSpec, mpnn_selection
from .Ui_MpnnViewer import Ui_MpnnViewer

KV_RE = re.compile(r"(?P<key>\w+)=(?P<value>[\w\.]+)")

def render_header_fn(record: SeqRecord, is_group_header: bool) -> str:

    if not is_group_header:
        return "/"
    
    attributes = dict(
        (m.group('key'), m.group('value'))
        for m in KV_RE.finditer(record.description)
    )

    sample = attributes.get('id')
    record_id = record.id.replace(",", "") if record.id is not None else "<unknown>"

    # if sample=X is not present, then this is the sequnece
    # of the structure
    if sample is None:
        return str(record_id)
    else:
        t = attributes.get('T')
        score = attributes.get('overall_confidence')
        return f"{record_id}, id={sample}, T={t}, score={score}"

def clean_id(id: Optional[str]) -> Optional[str]:

    if id is None:
        return None
    
    return id.replace(",", "").upper()

def get_selection(record: SeqRecord) -> Optional[List[StructureSelection]]:
    record_id = clean_id(record.id)

    if record_id is None:
        return None

    models: List[str] = cmd.get_names()
    model_ix = next(
        (
            ix
            for record in record_id.split(",")
            for ix, candidate in enumerate(models)
                if candidate.upper() == record.upper()
        ),
        None
    )

    if model_ix is None:
        return None
    
    model = models[model_ix]
    available_chains : List[str] = [
        chain.upper()
        for chain in cmd.get_chains(f"model {model} & polymer.protein & backbone")
    ]
    available_chains.sort()

    result = [
        StructureSelection(
            structure_name=model,
            chain_name=chain,
            segment_identifier=None
        )
        for chain in available_chains
    ]

    if len(result) > 0:
        return result
    else:
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
        self.__resv_map = spec.get_residues()
        self.__ui = Ui_MpnnViewer()
        self.__ui.setupUi(self)
        self.__processess : Optional[ExecutableProcessGroup] = None
        self.__ui.resultsFolderButton.clicked.connect(self.__on_results_folder_clicked)
        self.__ui.fastaViewer.set_render_header(render_header_fn)
        self.__ui.hideUndesignedCheckbox.stateChanged.connect(self.__on_hide_undesigned_checked)

        self.__init_widget()

    @pyqtSlot(int)
    def __on_hide_undesigned_checked(self, state: int):
        self.__update_visibility()

    @pyqtSlot()
    def __on_results_folder_clicked(self):
        QDesktopServices.openUrl(
            QUrl.fromLocalFile(self.__working_directory)
        )
    
    def __get_assigned_chains_jonsl_path(self) -> str:
        return path.join(self.__working_directory, "assigned.jsonl")
    
    def __get_fixed_positions_jonsl_path(self) -> str:
        return path.join(self.__working_directory, "fixed.jsonl")
    
    def __get_excluded_residues_jonsl_path(self) -> str:
        return path.join(self.__working_directory, "excluded.jsonl")
    
    def __get_tied_jonsl_path(self) -> str:
        return path.join(self.__working_directory, "tied.jsonl")
    
    def __get_model_location(self, model: str) -> str:
        return path.join(self.__working_directory, f"{model}.pdb")
    
    def __get_results_location(self) -> str:
        return path.join(self.__working_directory, "results")

    def __get_logfile_location(self) -> str:
        return path.join(self.__working_directory, "log.txt")
    
    def __get_model_locations_json(self, models: Sequence[str]) -> str:
        location = path.join(self.__working_directory, "models.json")

        with open(location, 'w') as location_stream:
            json.dump(
                {
                    self.__get_model_location(model): ''
                    for model in models
                },
                location_stream
            )

        return location

    def __save_models(self, models: Iterable[str]):

        for model in models:
            cmd.save(
                self.__get_model_location(model),
                mpnn_selection(model, include_ligands=self.__spec.mpnn_model.include_ligand)
            )

    def __models_to_pdb(self):
        return {
            model: self.__get_model_location(model)
            for model in self.__spec.get_models()
        }

    def __save_fixed_chains(self):
        spec = self.__spec
        with open(self.__get_assigned_chains_jonsl_path(), 'w') as jonsl:
            json.dump(spec.get_chains_jsonl(self.__models_to_pdb()), jonsl)

    def __save_fixed_positions(self):
        spec = self.__spec
        with open(self.__get_fixed_positions_jonsl_path(), 'w') as jsonl:
            json.dump(spec.get_positions_jsonl(self.__models_to_pdb()), jsonl)

    def __populate_fasta(self):

        seqs: List[SeqRecord] = []
        selections: Dict[int, SequenceStructures] = {}
        for fasta in Path(self.__get_results_location()).rglob("*.[Ff][aA]"):
            current = list(cast(Iterable[SeqRecord], SeqIO.parse(fasta, format='fasta')))

            if len(current) < 0:
                continue

            main = current[0]
            start_ix = len(seqs)
            seqs += current
            
            # Find out which of the currently loaded structures was used
            # as a template to generate the given sequences
            selection: Optional[List[Optional[StructureSelection]]] = cast(Any, get_selection(main))
            if selection is None:
                continue

            for i in range(start_ix, len(seqs)):
                selections[i] = selection

        self.__ui.fastaViewer.set_sequences(seqs, structure_mappings=selections)
        self.__update_visibility()

    def __update_visibility(self):
        hide_undesigned = self.__ui.hideUndesignedCheckbox.isChecked()
        sequences = self.__ui.fastaViewer.sequences()
        size = max(len(seq) for seq in sequences)

        if not hide_undesigned:
            mask = [
                True
                for _ in range(size)
            ]
        else:
            mask = [False for _ in range(size)]

            for space in self.__spec.edit_spaces:
                mappings = self.__resv_map[space.backbone]
                for resi in space.residues:
                    pos = mappings[resi].seq_pos
                    mask[pos] = True

        self.__ui.fastaViewer.set_fasta_position_mask(mask)

    @pyqtSlot(object)
    def __on_complete(self, result: ExecutableGroupResult):

        processes = self.__processess

        if processes is None:
            return

        processes.on_complete.disconnect(self.__on_complete)
        processes.close()
        commands = processes.get_commands()
        processes = None

        self.__ui.progressBar.hide()

        errors = result.get_errors()
        if errors is not None:
            show_error(self, "Protein MPNN Failed", errors)

        self.__populate_fasta()

        self.__write_log("\n".join(commands), result)

    def __write_log(self, command: str, result: ExecutableGroupResult) -> None:
        with open(self.__get_logfile_location(), 'w') as out_stream:
            out_stream.write(f"Command:\n{command}\n")
            out_stream.write(f"StdOut:\n {result.get_outputs()}\n")
            out_stream.write(f"Errors:\n{result.get_errors()}\n")

    def __get_excluded_args(self) -> List[str]:
        excluded = self.__spec.get_excluded_jonsl(self.__models_to_pdb())
        if len(excluded) < 1:
            return []
        
        excluded_file = self.__get_excluded_residues_jonsl_path()
        with open(excluded_file, 'w') as out_stream:
            json.dump(excluded, out_stream)

        return ["--omit_AA_per_residue_multi", excluded_file]
    
    def __get_tied_jonsl(self) -> Optional[str]:

        tied_jonsl = self.__spec.get_tied_jonsl(self.__models_to_pdb())

        if len(tied_jonsl) == 0:
            return None
        
        tied_file = self.__get_tied_jonsl_path()
        with open(tied_file, 'w') as jsonl:
            json.dump(tied_jonsl, jsonl)

        return tied_file

    def __run_mpnn(self):

        args = self.__spec.mpnn_args
        models = list(self.__spec.get_models())

        tied_jonsl = self.__get_tied_jonsl()
        if tied_jonsl is not None:
            tied_args = ["--symmetry_residues_multi", tied_jonsl]
        else:
            tied_args = []

        self.__processess = ExecutableProcessGroup.create([
            ExecutableProcess.create_process(
                self.__mpnn,
                [
                    #"--pdb_path", self.__get_model_location(model),
                    "--out_folder", self.__get_results_location(),
                    "--model_type", self.__spec.mpnn_model.model_name,
                    "--pdb_path_multi", self.__get_model_locations_json(models),
                    "--chains_to_design_multi", self.__get_assigned_chains_jonsl_path(),
                    #"--jsonl_path", self.__get_chains_jonsl_path(),
                    #"--chain_id_jsonl", self.__get_assigned_chains_jonsl_path(),
                    "--fixed_residues_multi", self.__get_fixed_positions_jonsl_path(),
                    "--batch_size", "1",
                    "--number_of_batches", str(self.__spec.num_seqs),
                    #"--num_seq_per_target", str(self.__spec.num_seqs),
                    #"--backbone_noise", str(args.backbone_noise),
                    "--temperature", str(args.sampling_temperature),
                    "--seed", str(args.random_seed),
                    "--verbose", "1",
                    "--save_stats", "1"
                ] \
                + self.__get_excluded_args() \
                + tied_args
            )
        ])

        self.__processess.on_complete.connect(self.__on_complete)
        self.__processess.start()

    def __init_widget(self):
        self.__save_models(self.__spec.get_models())
        self.__save_fixed_chains()
        self.__save_fixed_positions()
        self.__run_mpnn()




