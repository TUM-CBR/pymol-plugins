from io import StringIO, UnsupportedOperation
import subprocess
from typing import List, Optional, TextIO

from . import streams

def only_if_file(stream : TextIO) -> Optional[TextIO]:

    try:
        stream.fileno()
        return stream
    except UnsupportedOperation:
        return None

def simple_execute(
    command : List[str],
    input : streams.InputStream = "",
    output : Optional[TextIO] = None,
    cwd : Optional[str] = None
    ) -> None:

    with streams.get_input_stream(input) as input_io \
        , StringIO() as errors \
        , subprocess.Popen(
            command,
            text=True,
            cwd = cwd,
            stdin = only_if_file(input_io.stream) or subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        ) as process:

        if process.stdin:
            for text in input_io.stream:
                process.stdin.write(text)
            process.stdin.close()

        if process.stdout:
            for text in process.stdout:
                 if output:
                    output.write(text)

        if process.stderr:
            for text in process.stderr:
                errors.write(text)

        if process.wait() != 0:
            errors.seek(0)
            raise IOError(errors.read())