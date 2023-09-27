from io import StringIO
from os import path
import os
import subprocess
import tempfile
from typing import List, NamedTuple, TextIO

from ..support.msa import Msa

ACPSICOV_NOT_FOUND="""Only windows binaries are currently provided with 'cbr-tools'. The feature you
are trying to use requires 'ACPSICOV'. Plese obtain a copy at https://github.com/etaijacob/ACPSICOV
and ensure it is in your PATH.
"""

def acpsicov_bin():

    if os.name == 'nt':
        return os.path.join(
            os.path.dirname(__file__),
            "resources",
            "amc.exe"
        )

    else:
        import shutil
        if shutil.which("amc"):
            return "amc"
        else:
            raise Exception(ACPSICOV_NOT_FOUND)

class AcpsicovEntry(NamedTuple):
    position1 : int
    position2 : int
    confidence : float

    def write(self, stream : TextIO):
        stream.write(f"{self.position1} {self.position2} {self.confidence}\n")

class AcpsicovResult(NamedTuple):
    entries : List[AcpsicovEntry]

    def write(self, stream: TextIO):

        stream.write("i j psicov.aa\n")

        for entry in self.entries:
            entry.write(stream)

def acpsicov_load(stream : TextIO) -> AcpsicovResult:

    return AcpsicovResult(
        entries = [
            AcpsicovEntry(position1=int(items[0]), position2=int(items[1]), confidence=float(items[2]))
            for line in stream.readlines()[1:]
            for items in [line.split()]
        ]
    )

def acpsicov(msa : Msa) -> AcpsicovResult:

    with tempfile.TemporaryDirectory() as temp_dir:
        seqs_file_name = path.join(temp_dir, "sequences.txt")
        prev_length = None

        with open(seqs_file_name, 'w') as seqs_file:

            for sequence in msa.values():

                if prev_length is not None and prev_length != len(sequence):
                    raise ValueError(
                        f"Dodgey MSA. The length of all sequences must be identical."
                    )
                prev_length = len(sequence)

                seqs_file.write(f"{sequence}\n")

        process = subprocess.Popen(
            [acpsicov_bin(), "-o", "-p", "-d", "0.03", seqs_file_name],
            text=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        with \
            process, \
            StringIO() as out_string, \
            StringIO() as out_err:

            assert process.stdout, "Could not open standard output of amc"
            for text in process.stdout:
                out_string.write(text)

            assert process.stderr, "Could not open standard error of amc"
            for text in process.stderr:
                out_err.write(text)

            if process.wait() != 0:
                out_err.seek(0)
                error_str = out_err.read()
                raise Exception(error_str)

            out_string.seek(0)
            return acpsicov_load(out_string)