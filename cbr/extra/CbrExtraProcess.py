from io import StringIO
import json
from PyQt5.QtCore import QIODevice, QObject, QProcess, pyqtSignal, pyqtSlot
import os
import subprocess
from typing import Any, Callable, Dict, List, NamedTuple, Optional, Sequence, TextIO, TypeVar

NOT_FOUND_ERROR="""Only windows binaries are currently provided with 'cbr-tools'. The feature you
are trying to use requires 'cbr-tools-extra'. Plese obtain a copy at https://github.com/TUM-CBR/cbr-tools-extra
and ensure it is in your PATH.
"""

def cbrtools_bin() -> str:

    if os.name == 'nt':
        return os.path.join(
            os.path.dirname(__file__),
            "resources",
            "cbrtools",
            "cbrtools.exe"
        )

    else:
        import shutil
        cbrtools_cmd = shutil.which("cbrtools")
        if cbrtools_cmd is not None:
            return cbrtools_cmd
        else:
            raise Exception(NOT_FOUND_ERROR)

T = TypeVar('T')

class CbrToolsProcessException(Exception):
    pass

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

class Logger:

    def __init__(self, log_location: Optional[str]):
        self.__session_log = None if log_location is None else open(log_location, 'w')

    def log_message(self, message: str):

        if self.__session_log is None:
            return
        
        newline = "\n" + "-"*100 + "\n"
        self.__session_log.write(f"{message}{newline}")
        self.__session_log.flush()

    def close(self):
        if self.__session_log is None:
            return

        self.__session_log.close()

class CbrExtraProcess(QProcess):

    message_signal = pyqtSignal(object)
    error_signal = pyqtSignal(object)

    def __init__(self, *args : Any, **kwargs: Any):
        super().__init__(*args, **kwargs)

        self.setProgram(cbrtools_bin())
        self.readyReadStandardOutput.connect(self.__on_data_ready)
        self.finished.connect(self.__on_process_finish)

        self.__logger_instance : Optional[Logger] = Logger(os.environ.get('LOG_CBR_EXTRA'))

    @property
    def __logger(self) -> Logger:
        result = self.__logger_instance

        assert result is not None, "Logger requested after close"
        return result

    def __close_logger(self):
        logger = self.__logger_instance

        if logger is None:
            return
        
        self.__logger_instance = None
        logger.close()

    @pyqtSlot(int, QProcess.ExitStatus)
    def __on_process_finish(self, exit_code: int, exit_status: QProcess.ExitStatus):

        self.__close_logger()

        if exit_code == 0:
            return
        self.setCurrentReadChannel(QProcess.ProcessChannel.StandardError)
        result = decode_bytes(self.readAll())
        self.error_signal.emit(Exception(result))

    def setProcessChannelMode(self, mode: int) -> None:
        raise Exception("The read channel should not be changed by the user")

    def setReadChannel(self, channel: 'QProcess.ProcessChannel') -> None:
        raise Exception("The read channel should not be changed by the user")

    def run_cbr_process(self, args : List[str]):
        self.setArguments(args)
        self.start(QIODevice.ReadWrite | QIODevice.Text) # type: ignore

    def write_json_dict(self, value: Dict[Any, Any]):
        message = f"{json.dumps(value)}\n"
        self.__logger.log_message(message)
        self.write(message.encode('utf-8'))

    def close(self):
        self.write_json_dict({
            'uid': -1,
            'entity_type': 'stop'
        })

        self.__close_logger()
        if not super().waitForFinished(msecs=1000):            
            super().close()

    @pyqtSlot()
    def __on_data_ready(self):

        assert self.readChannel() == QProcess.ProcessChannel.StandardOutput, "Read channel should be stdin"

        while self.canReadLine():
            line_bytes : Any = self.readLine()
            line_str = decode_bytes(line_bytes)
            self.__logger.log_message(line_str)

            try:
                value = json.loads(line_str)
            except json.JSONDecodeError:
                continue

            self.message_signal.emit(value)

def run_cbr_tools_interactive(
    args : List[str],
    parent : Optional[QObject] = None
) -> CbrExtraProcess:

    process = CbrExtraProcess(parent=parent)
    process.run_cbr_process(args)

    return process

def run_cbr_tools(
    args : List[str],
    in_data : Optional[TextIO] = None,
    out_handler : Optional[Callable[[TextIO], T]] = None
) -> Optional[T]:

    process = subprocess.Popen(
        [cbrtools_bin()] + args,
        text = True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    with process, \
        StringIO() as error_stream, \
        StringIO() as out_stream:

        assert process.stdin, "Could not open the standard output of 'cbrtools'"

        if in_data:
            for data in in_data:
                process.stdin.write(data)
            process.stdin.close()

        assert process.stdout, "Could not open the standard output of 'cbrtools'"
        for data in process.stdout:
            out_stream.write(data)

        assert process.stderr, "Could not open the standard error of 'cbrtools'"
        for text in process.stderr:
            error_stream.write(text)

        if process.wait() != 0:
            error_stream.seek(0)
            error_str = error_stream.read()
            raise CbrToolsProcessException(error_str)

        if out_handler is not None:
            out_stream.seek(0)
            return out_handler(out_stream)

    return None

class CommandResult(NamedTuple):
    std_out: StringIO
    std_error: StringIO
    exit_code: int
    command_id: Any = None

    def __enter__(self, *args: Any, **kwargs: Any) -> 'CommandResult':
        self.std_out.__enter__(*args, **kwargs)
        self.std_error.__enter__(*args, **kwargs)
        return self
    
    def __exit__(self, *args: Any, **kwargs: Any):
        self.std_out.__exit__(*args, **kwargs)
        self.std_error.__exit__(*args, **kwargs)

class CommandInput(NamedTuple):
    command_id: Any
    args: Sequence[str]
    input: Optional[TextIO]
    env: Optional[Dict[str, str]] = None

class CBRCommandRunner(QObject):

    command_done_signal = pyqtSignal(object)
    command_run_signal = pyqtSignal(object)

    def __init__(self, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)

        self.command_run_signal.connect(self.__on_run_process)

    @pyqtSlot(object)
    def __on_run_process(self, cmd: CommandInput):

        std_err = StringIO()
        std_out = StringIO()

        if os.name == "nt":
            flags = subprocess.CREATE_NO_WINDOW
        else:
            flags = 0

        args: List[str] = [cbrtools_bin()] + list(cmd.args)
        print(args)

        cbr_process = subprocess.Popen(
            args,
            text = True,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
            creationflags = flags,
            env = cmd.env
        )

        std_in = cmd.input if cmd.input is not None else StringIO()

        with cbr_process, std_in:

            assert cbr_process.stdin, "Standard input is expected to be open."
            for text in  std_in:
                cbr_process.stdin.write(text)
            cbr_process.stdin.close()

            assert cbr_process.stdout, "Standard output is expected to be open"
            for text in cbr_process.stdout:
                std_out.write(text)

            assert cbr_process.stderr, "Standard error is expected to be open"
            for text in cbr_process.stderr:
                std_err.write(text)

            cbr_process.wait()

            result = cbr_process.returncode

        std_out.seek(0)
        std_err.seek(0)
        self.command_done_signal.emit(
            CommandResult(
                std_out,
                std_err,
                result,
                command_id=cmd.command_id
            )
        )

    def __default_env__(self) -> Optional[Dict[str, str]]:
        return None

    def __create_env(self, env: Optional[Dict[str, str]]) -> Dict[str, str]:

        return {
            k:v
            for current_env in [self.__default_env__(), env]
                if current_env is not None
            for k,v in current_env.items()
        }

    def run_command(
        self,
        args: Sequence[str],
        input: Optional[TextIO] = None,
        command_id: Any = None,
        env: Optional[Dict[str, str]] = None
    ):
        
        self.command_run_signal.emit(
            CommandInput(
                args = args,
                input = input,
                command_id = command_id,
                env=self.__create_env(env)
            )
        )