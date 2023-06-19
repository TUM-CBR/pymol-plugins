from io import StringIO, TextIOBase
import os
import subprocess
from typing import Callable, Iterable, NewType, TextIO, Tuple, TypeVar

from . import fasta

def find_clustal():

    if os.name == 'nt':
        return path.join(path.dirname(__file__), "resources", "clustalo.exe")

    return "clustalo"

class ClustalResult(object):
    pass

class WrapIO(object):

    def __init__(
        self,
        stream : 'TextIO | None' = None,
        open_stream : 'Callable[[], TextIO]' = StringIO,
        init : 'None | Callable[[TextIO], None]' = None):
        self.__stream = stream or open_stream()
        self.__init = init
        self.__percolate = stream and True

    @property
    def stream(self) -> TextIO:
        return self.__stream

    def __enter__(self, *args, **kwargs):

        if self.__percolate:
            self.__stream.__enter__(*args, **kwargs)

        if self.__init:
            self.__init(self.__stream)

        return self

    def __exit__(self, *args, **kwargs):

        if self.__percolate:
            self.__stream.__exit__(*args, **kwargs)

MsaInput = TypeVar('MsaInput', str, Iterable[Tuple[str, str]], TextIO)

MsaOutput = TypeVar('MsaOutput', str, TextIO)

class Clustal(object):

    def __init__(
        self,
        clustal_executable = None
    ):

        self.__clustal_executable = clustal_executable or find_clustal()

    @staticmethod
    def __get_msa_input(in_msa : MsaInput) -> WrapIO:

        if(isinstance(in_msa, TextIO)):
            return WrapIO(stream = in_msa)

        def init(result : TextIO):
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

        if(isinstance(out_msa, TextIO)):
            return WrapIO(stream = out_msa)
        elif(isinstance(out_msa, str)):
            return WrapIO(open_stream=lambda: open(out_msa, 'w'))
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

            process.wait()

        return ClustalResult()

    def run_msa_items(self, items : MsaInput) -> 'Iterable[Tuple[str, str] | Exception]':

        with StringIO() as result:
            self.run_msa(items, result)
            return fasta.parse_fasta_stream(result)