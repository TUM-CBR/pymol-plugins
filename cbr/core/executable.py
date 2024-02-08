from enum import Enum
from PyQt5.QtCore import QObject, QProcess, QTimer, pyqtSignal, pyqtSlot
from typing import Any, Iterable, List, NamedTuple, Optional

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

class ExecutableResult(NamedTuple):
    process :'ExecutableProcess'
    error : Optional[str]

    @property
    def is_successful(self) -> bool:
        return self.error is not None

class ExecutableProcess(QProcess):
    """
    Base class to run external processess and monitor their
    execution. Offers simple signals to get the result
    of running a program.
    """

    on_completed = pyqtSignal(object)

    def __init__(
        self,
        executable: Executable,
        args: List[str]
    ):
        super().__init__()
        self.setProgram(executable.location)
        self.setArguments(args)
        self.finished.connect(self.__on_process_finish)
        self.__executable_result : Optional[ExecutableResult] = None

    def get_result(self) -> Optional[ExecutableResult]:

        result = self.__executable_result
        if result is not None:
            return result
        
        state = self.state()

        if state != QProcess.ProcessState.NotRunning:
            return None

        exit_code = self.exitCode()

        if exit_code == 0:
            result = ExecutableResult(
                process = self,
                error = None
            )
        else:
            result = ExecutableResult(
                process = self,
                error = decode_bytes(self.readAllStandardError())
            )

        self.__executable_result = result

        return result

    @pyqtSlot(int, QProcess.ExitStatus)
    def __on_process_finish(self, exit_code: int, exit_status: QProcess.ExitStatus):

        result = self.get_result()
        self.on_completed.emit(result)

    @staticmethod
    def create_process(executable: Executable, args: List[str]) -> 'ExecutableProcess':
        process = ExecutableProcess(executable, args)
        return process
    
class ExecutableGroupResult(NamedTuple):
    results : List[ExecutableResult]

    def get_errors(self) -> Optional[str]:
        errors = [
            result.error
            for result in self.results if result.error is not None
        ]

        if len(errors) == 0:
            return None
        
        return "\n".join(errors)

class ExecutableProcessGroup(QObject):
    """
    Class used for managing multiple processess. It assists in keeping track
    of subscriptions to multiple processess and issuing signals based on the
    results of multiple processess.

    Attributes:
    on_complete : Signal[ExecutableGroupResult]
    """

    on_complete = pyqtSignal(object)

    def __init__(
        self,
        processess: Iterable[ExecutableProcess]
    ):
        super().__init__()
        self.__processess = list(processess)
        self.__ping_timer = QTimer()
        self.__ping_timer.timeout.connect(self.__on_timeout)

        self.__connections = [
            process.on_completed.connect(self.__on_complete)
            for process in processess
        ]

    def __conclude(self):

        self.__ping_timer.stop()

        def assert_result(result: Optional[ExecutableResult]) -> ExecutableResult:

            assert result is not None, "This method should not be called unless all processess have concluded"
            return result

        self.on_complete.emit(
            ExecutableGroupResult(
                results=[
                    assert_result(process.get_result())
                    for process in self.__processess
                ]
            )
        )

    def close(self):

        for process in self.__processess:
            process.close()

        for connection in self.__connections:
            self.disconnect(connection)

    @pyqtSlot()
    def __on_complete(self, result: ExecutableResult):
        self.__ping(result)

    @pyqtSlot()
    def __on_timeout(self):
        self.__ping(None)

    def __ping(self, _result: Optional[ExecutableResult]):

        not_running = (
            process.state() == QProcess.ProcessState.NotRunning
            for process in self.__processess
        )

        if all(not_running):
            self.__conclude()

    def start(self):

        for process in self.__processess:
            process.start()

        ping = self.__ping_timer
        ping.setInterval(1000)
        ping.start()

    @staticmethod
    def create(members: Iterable[ExecutableProcess]) -> 'ExecutableProcessGroup':
        return ExecutableProcessGroup(members)