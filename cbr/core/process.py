from io import StringIO
import subprocess
from typing import List, Optional, TextIO

from . import streams

def simple_execute(
    command : List[str],
    input : streams.InputStream = "",
    output : Optional[TextIO] = None
    ) -> None:

    process = subprocess.Popen(
        command,
        text=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    with process \
        , streams.get_input_stream(input) as input_io \
        , StringIO() as errors:

        if process.stdin:
            for text in input_io.stream:
                process.stdin.write(text)
                process.stdin.close()
        else:
            raise Exception("Failed to open standard input for clustal")

        if process.stdout:
            for text in process.stdout:
                 if output:
                    output.write(text)

        if process.stderr:
            for text in process.stderr:
                errors.write(text)

        if process.wait() != 0:
            raise IOError(errors.read())