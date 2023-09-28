from io import StringIO
import json
import subprocess
from typing import Callable, List, Optional, TextIO, TypeVar

from ..extra.main import cbrtools_bin
from .data import DesignPrimersResults

T = TypeVar('T')

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
            raise Exception(error_str)

        if out_handler is not None:
            out_stream.seek(0)
            return out_handler(out_stream)

    return None

def design_primers(
    start : int,
    codon_count : int,
    min_length : int,
    max_length : int,
    organism : str,
    sequence : str,
    result_db : str
):
    with StringIO() as in_data:
        json.dump(
            {
                "start": start,
                "codon_count": codon_count,
                "min_length": min_length,
                "max_length": max_length,
                "organism": organism,
                "sequence": sequence
            },
            in_data
        )
        in_data.seek(0)
        run_cbr_tools(
            ["primers", "design", result_db],
            in_data=in_data
        )

def query_best_primers(
    db : str,
    tm : float,
    wTm : float,
    pTm : float,
    wPTm : float,
    wTmDelta : float
) -> DesignPrimersResults:

    result = run_cbr_tools(
        [
            "primers",
            "query",
            f"--tm_primers={pTm},{wPTm}",
            f"--tm_total={tm},{wTm}",
            f"--tm_delta={wTmDelta}",
            db
        ],
        out_handler=lambda out_stream: DesignPrimersResults.from_json(json.load(out_stream))
    )

    assert result is not None, "Bug in the code, 'run_cbr_tools' returned None with a handler."
    return result