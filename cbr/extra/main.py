from io import StringIO
import json
from PyQt5.QtCore import QIODevice, QProcess, pyqtSignal, pyqtSlot
import os
import subprocess
from typing import Any, Callable, List, Optional, TextIO, TypeVar

NOT_FOUND_ERROR="""Only windows binaries are currently provided with 'cbr-tools'. The feature you
are trying to use requires 'cbr-tools-extra'. Plese obtain a copy at https://github.com/TUM-CBR/cbr-tools-extra
and ensure it is in your PATH.
"""

def cbrtools_bin():

    if os.name == 'nt':
        return os.path.join(
            os.path.dirname(__file__),
            "resources",
            "cbrtools",
            "cbrtools.exe"
        )

    else:
        import shutil
        if shutil.which("cbrtools"):
            return "cbrtools"
        else:
            raise Exception(NOT_FOUND_ERROR)

T = TypeVar('T')

class CbrToolsProcessException(Exception):
    pass

def decode_bytes(bs):
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

class CbrExtraProcess(QProcess):

    message_signal = pyqtSignal(object)
    error_signal = pyqtSignal(object)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.setProgram(cbrtools_bin())
        self.readyReadStandardOutput.connect(self.__on_data_ready)
        self.finished.connect(self.__on_process_finish)

        self.__log_values = os.environ.get('LOG_CBR_EXTRA')

    @pyqtSlot(int, QProcess.ExitStatus)
    def __on_process_finish(self, exit_code, exit_status):

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
        self.start(QIODevice.ReadWrite | QIODevice.Text)

    def write_json_dict(self, value: dict):
        bs = f"{json.dumps(value)}\n".encode('utf-8')
        self.write(bs)

    def close(self):
        self.write_json_dict({
            'uid': -1,
            'entity_type': 'stop'
        })
        super().close()

    @pyqtSlot()
    def __on_data_ready(self):

        assert self.readChannel() == QProcess.ProcessChannel.StandardOutput, "Read channel should be stdin"

        while self.canReadLine():
            line_bytes : Any = self.readLine()

            if self.__log_values:
                print(line_bytes)

            try:
                value = json.loads(decode_bytes(line_bytes))
            except json.JSONDecodeError:
                continue

            self.message_signal.emit(value)

def run_cbr_tools_interactive(
    args : List[str]
) -> CbrExtraProcess:

    process = CbrExtraProcess()
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