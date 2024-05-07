from PyQt5.QtCore import pyqtSignal
from typing import NamedTuple, Optional, Sequence

from ..extra.CbrExtraProcess import CBRCommandRunner

class RunScanResult(NamedTuple):
    error: Optional[str] = None

class SequenceCommandRunner(CBRCommandRunner):

    scan_done_signal = pyqtSignal()

    def run_scan(
        self,
        db_file: str,
        scan_folders: Optional[Sequence[str]] = None
    ):
    