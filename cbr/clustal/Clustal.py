from io import StringIO, TextIOBase
import os
import subprocess
from typing import cast, Dict, Iterable, Tuple, TypeVar

from ..core.Context import Context
from ..core.WrapIO import WrapIO
from . import msa

def find_clustal():

    if os.name == 'nt':
        return os.path.join(os.path.dirname(__file__), "resources", "clustalo.exe")

    return "clustalo"

class ClustalResult(object):
    pass

MsaInput = TypeVar('MsaInput', str, Iterable[Tuple[str, str]], TextIOBase)

MsaOutput = TypeVar('MsaOutput', str, TextIOBase)

def get_clustal() -> 'Clustal':
    return Clustal()

def get_clustal_from_context(context : Context):
    return Clustal()

class Clustal(object):

    def __init__(
        self,
        clustal_executable = None
    ):

        self.__clustal_executable = clustal_executable or find_clustal()

    @staticmethod
    def __get_msa_input(in_msa : MsaInput) -> WrapIO:

        if(isinstance(in_msa, TextIOBase)):
            return WrapIO(stream = in_msa)

        def init(result : TextIOBase):
            if isinstance(in_msa, str):
                result.write(in_msa)

            elif isinstance(in_msa, Iterable):
                for (name, sequence) in in_msa:
                    result.write('>%s\n%s\n\n' % (name, sequence))
            else:
                raise ValueError('The input to clustal is not valid')

        return WrapIO(init = init)

    @staticmethod
    def __get_msa_output(out_msa : MsaOutput) -> WrapIO:

        if(isinstance(out_msa, TextIOBase)):
            return WrapIO(stream = cast(TextIOBase, out_msa))
        elif(isinstance(out_msa, str)):
            return WrapIO(open_stream = lambda: open(out_msa, 'w'))
        else:
            raise ValueError('The output to clustal is not valid')

    def run_msa(
        self,
        in_msa : MsaInput,
        out_result : MsaOutput
        ) -> ClustalResult:

        process = subprocess.Popen(
            [self.__clustal_executable, '-i', '-', '--force', '--outfmt=clustal'],
            text=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        with process \
            , Clustal.__get_msa_input(in_msa) as in_msa_stream \
            , Clustal.__get_msa_output(out_result) as out_result_stream:

            if process.stdin:
                for text in in_msa_stream.stream:
                    process.stdin.write(text)
                process.stdin.close()
            else:
                raise Exception("Failed to open standard input for clustal")
            
            if process.stdout:
                for text in process.stdout:
                    out_result_stream.stream.write(text)
            else:
                raise Exception("Failed to open the standard output of clustal")

            if process.wait() != 0:
                raise ValueError("Invalid input provided to clustal.")
            

        return ClustalResult()

    def run_msa_items(self, items : MsaInput) -> Dict[str, str]:

        with StringIO() as result:
            self.run_msa(items, result)
            result.seek(0)
            return msa.parse_alignments(result)
