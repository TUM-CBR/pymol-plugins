import tempfile
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO, TextIOBase
import os
import subprocess
from typing import Any, cast, Dict, Iterable, List, NamedTuple, Optional, Tuple, TypeVar

from ..core.Context import Context
from ..core.WrapIO import WrapIO
from ..core import assertions
from . import msa

def find_clustal():

    if os.name == 'nt':
        return os.path.join(os.path.dirname(__file__), "resources", "clustalo.exe")

    return "clustalo"

class ClustalResult(object):
    pass

class MsaFromFragments(NamedTuple):
    alignment : MultipleSeqAlignment
    fragments_positions : List[int]

MsaInput = TypeVar('MsaInput', str, Iterable[Tuple[str, str]], TextIOBase)

MsaOutput = TypeVar('MsaOutput', str, TextIOBase)

def get_clustal() -> 'Clustal':
    return Clustal()

def get_clustal_from_context(context : Context):
    return Clustal()

class Clustal(object):

    def __init__(
        self,
        clustal_executable : Optional[str] = None
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

        out_dir = tempfile.TemporaryDirectory()
        out_file = os.path.join(out_dir.name, "out.clustal")

        if os.name == "nt":
            flags = subprocess.CREATE_NO_WINDOW
        else:
            flags = 0

        process = subprocess.Popen(
            [
                self.__clustal_executable,
                '-i',
                '-',
                '--force',
                '--outfmt=clustal',
                '-o',
                out_file
            ],
            text=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            creationflags=flags
        )

        with process \
            , out_dir \
            , Clustal.__get_msa_input(in_msa) as in_msa_stream \
            , StringIO() as error_stream \
            , Clustal.__get_msa_output(out_result) as out_result_stream:

            if process.stdin:
                for text in in_msa_stream.stream:
                    process.stdin.write(text)
                process.stdin.close()
            else:
                raise Exception("Failed to open standard input for clustal")
            
            process.wait()

            with open(out_file, 'r') as fs:
                for text in fs.readlines():
                    out_result_stream.stream.write(text)

            if process.stderr is not None:
                for text in process.stderr:
                    error_stream.write(text)
            else:
                error_stream.write("Invalid input provided to clustal.")

            if process.wait() != 0:
                error_stream.seek(0)
                raise ValueError(error_stream.read())
            

        return ClustalResult()
    
    def run_msa_seqs(self, seqs: List[SeqRecord]) -> MultipleSeqAlignment:
        with StringIO() as seqs_in, \
            StringIO() as seqs_out:
            SeqIO.write(
                seqs,
                seqs_in,
                format='fasta'
            )
            seqs_in.seek(0)
            self.run_msa(seqs_in, seqs_out)
            seqs_out.seek(0)
            seq_it: Any = AlignIO.parse(
                seqs_out,
                format='clustal'
            )
            return next(seq_it)

    def run_msa_items(self, items : MsaInput) -> Dict[str, str]:

        with StringIO() as result:
            self.run_msa(items, result)
            result.seek(0)
            return msa.parse_alignments(result)
        
    def run_msa_fragments(
        self,
        fragments: Dict[str, List[str]],
        result_order : Optional[List[str]] = None
    ) -> MsaFromFragments:
        
        if result_order is None:
            result_order = list(fragments.keys())

        fragment_entries = [
            (parts[0], " ".join(parts[1:]), fragments[item])
            for item in result_order
            for parts in [item.split()]
        ]

        count = assertions.assert_same_length(*(list(values) for values in fragments.values()))
        result : Dict[str, str] = dict(
            (id, "")
            for id,_,_ in fragment_entries
        )

        position = 0
        positions : List[int] = []
        for i in range(0, count):
            input = (
                (id, seq)
                for id,_,seqs in fragment_entries
                for seq in [seqs[i].strip()] if len(seq) > 0
            )
            partial_msa = self.run_msa_items(input)

            # Append the MSA results for the fragment to the
            # complete result
            for name,seq in partial_msa.items():
                result[name] += seq

            msa_length = assertions.assert_same_length(*partial_msa.values())
            position += msa_length

            # The position of the last fragment is the end of the sequence
            # this position is not needed
            if i + 1 < count:
                positions.append(position)

            # Add gaps for the sequences that were not included in the
            # alignment (ie. empty strings provided)
            for id,_,_ in fragment_entries:
                if id not in partial_msa:
                    result[id] += "-"*msa_length

        records = [
            SeqRecord(
                seq = Seq(result[id]),
                id = id,
                name = name
            )
            for id,name,_ in fragment_entries
        ]

        return MsaFromFragments(
            alignment = MultipleSeqAlignment(records),
            fragments_positions = positions
        )