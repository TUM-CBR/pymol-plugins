from Bio.Align import MultipleSeqAlignment
from os import path
from typing import Callable, Dict, TextIO

def fasta_writer(stream : TextIO, msa: MultipleSeqAlignment):
    stream.write(msa.__format__("fasta"))

def clustal_writer(stream : TextIO, msa: MultipleSeqAlignment):
    stream.write(msa.__format__("clustal"))

MSA_EXTENSIONS_WRITERS : Dict[str, Callable[[TextIO, MultipleSeqAlignment], None]] = {
    "fasta": fasta_writer,
    "fa": fasta_writer,
    "clustal": clustal_writer
}

def save_msa(location: str, msa : MultipleSeqAlignment):

    _, ext = path.splitext(location)

    writer = MSA_EXTENSIONS_WRITERS.get(ext.replace(".", ""))
    if writer is None:
        raise ValueError("The file {location} does not have a known msa format.")

    with open(location, 'w') as stream:
        writer(stream, msa)