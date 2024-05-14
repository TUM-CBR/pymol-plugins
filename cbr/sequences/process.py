from io import StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from enum import Enum
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot
from typing import Dict, Optional, Sequence
from warnings import warn

from ..blast.blast import Blast
from ..core import namedtuple
from ..extra.CbrExtraProcess import CBRCommandRunner, CommandResult

from .data import *

class RunScanResult(NamedTuple):
    error: Optional[str] = None

    @classmethod
    def from_command_result(cls, result: CommandResult) -> 'RunScanResult':

        if result.exit_code == 0:
            return RunScanResult()
        else:
            return RunScanResult(
                error = result.std_error.read()
            )
        
class RunSearchResult(NamedTuple):
    results: Optional[QueryResults] = None
    error: Optional[str] = None

    @classmethod
    def from_command_result(cls, result: CommandResult) -> 'RunSearchResult':

        if result.exit_code == 0:
            return RunSearchResult(
                results = namedtuple.load(QueryResults, result.std_out)
            )
        else:
            return RunSearchResult(error=result.std_error.read())

class RunErrorsResult(NamedTuple):
    results: Optional[ErrorResults] = None
    error: Optional[str] = None

    @classmethod
    def from_command_result(cls, result: CommandResult) -> 'RunErrorsResult':

        if result.exit_code == 0:
            return RunErrorsResult(
                results = namedtuple.load(ErrorResults, result.std_out)
            )
        else:
            return RunErrorsResult(error=result.std_error.read())

class CommandId(Enum):
    Scan = 0
    Query = 1
    Errors = 2

class SequenceCommandRunner(CBRCommandRunner):

    scan_done_signal = pyqtSignal(object)
    query_done_signal = pyqtSignal(object)
    errors_done_signal = pyqtSignal(object)

    def __init__(
        self,
        blast: Blast,
        parent: Optional[QObject] = None
    ) -> None:
        super().__init__(parent)

        self.__blast = blast
        self.command_done_signal.connect(self.__on_command_done)

    @pyqtSlot(object)
    def __on_command_done(self, result: CommandResult):

        with result:

            if result.command_id == CommandId.Scan:
                scan_result = RunScanResult.from_command_result(result)
                self.scan_done_signal.emit(scan_result)
            elif result.command_id == CommandId.Query:
                self.query_done_signal.emit(
                    RunSearchResult.from_command_result(result)
                )
            elif result.command_id == CommandId.Errors:
                self.errors_done_signal.emit(
                    RunErrorsResult.from_command_result(result)
                )
            else:
                warn(
                    f"The 'SequenceCommandRunner' got an unexpected command '{result.command_id}'."
                )


    def __default_env__(self) -> Optional[Dict[str, str]]:
        return {
            'CBR_MAKEBLAST_DB': self.__blast.makeblastdb,
            'CBR_TBLASTN': self.__blast.tblastn
        }

    def run_scan(
        self,
        db_file: str,
        scan_folders: Optional[Sequence[str]] = None
    ):
        
        if scan_folders is not None \
            and len(scan_folders) > 0:
            scan_folders_args = ["--"] + list(scan_folders)
        else:
            scan_folders_args = []


        self.run_command(
            [
                "sequences", "scan",
                "--db-file", db_file,
            ] + scan_folders_args,
            command_id = CommandId.Scan
        )

    def run_query(
        self,
        db_file: str,
        sequences: Sequence[SeqRecord]
    ):
        """
        Query the DNA database using cbrextra tools. This is a tblastx query
        meaning that a peptide fasta file is provided which will be used to
        search for genes coding the sequences in said files.
        """
        
        if len(sequences) < 1:
            raise ValueError("At least one sequence must be provided")
        
        input_seqs = StringIO()
        SeqIO.write(sequences, input_seqs, format='fasta')
        input_seqs.seek(0)

        self.run_command(
            [
                "sequences", "tblastn",
                "--db-file", db_file
            ],
            input = input_seqs,
            command_id = CommandId.Query
        )

    def run_query_errors(
        self,
        db_file: str
    ):
        self.run_command(
            [
                "sequences", "errors",
                "--db-file", db_file
            ],
            command_id = CommandId.Errors
        )