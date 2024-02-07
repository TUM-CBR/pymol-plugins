from enum import Enum
from PyQt5.QtCore import QIODevice, QProcess, pyqtSignal, pyqtSlot
from typing import Any, cast, List, NamedTuple

class ExecutableType(Enum):
    Binary = 0

class Executable(NamedTuple):
    executable_type : ExecutableType
    location : str

class KnownExecutables(Enum):
    ProteinMPNN = "ProteinMPNN"

def decode_bytes(bs: Any):
    try:
        # The typing claims that line_bytes should already
        # be bytes but experience shows that it is a
        # QByteArray, which should never be exposed to
        # python. LOoks like a PyQt5 bug, so we try both
        # and defend ourselves
        # Remember that with Python, types usually lie
        return bs.data().decode("utf-8")
    except AttributeError:
        return bs.decode("utf-8")

class ExecutableProcess(QProcess):
    """
    Base class to run external processess and monitor their
    execution. Offers simple signals to get the result
    of running a program.
    """

    success_signal = pyqtSignal()
    error_signal = pyqtSignal(object)

    def __init__(
        self,
        executable: Executable,
        args: List[str]
    ):
        super().__init__()
        self.setProgram(executable.location)
        self.setArguments(args)
        self.finished.connect(self.__on_process_finish)

    @pyqtSlot(int, QProcess.ExitStatus)
    def __on_process_finish(self, exit_code: int, exit_status: QProcess.ExitStatus):

        if exit_code == 0:
            self.success_signal.emit()

        self.setCurrentReadChannel(QProcess.ProcessChannel.StandardError)
        result = decode_bytes(self.readAll())
        self.error_signal.emit(Exception(result))

    @staticmethod
    def create_process(executable: Executable, args: List[str]) -> 'ExecutableProcess':
        process = ExecutableProcess(executable, args)
        return process