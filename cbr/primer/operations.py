from io import StringIO
import json
import subprocess

from ..extra.main import cbrtools_bin
from .data import DesignPrimersResults

def query_best_primers(
    db : str,
    tm : float,
    wTm : float,
    pTm : float,
    wPTm : float,
    wTmDelta : float
) -> DesignPrimersResults:

    process = subprocess.Popen(
        [
            cbrtools_bin(),
            "primers",
            "query",
            f"--tm_primers={pTm},{wPTm}",
            f"--tm_total={tm},{wTm}",
            f"--tm_delta={wTmDelta}",
            db
        ],
        text=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    with process, \
        StringIO() as out_string, \
        StringIO() as out_err:

        assert process.stdout, "Could not open the standard output of 'cbrtools'"
        for text in process.stdout:
            out_string.write(text)

        assert process.stderr, "Could not open the standard error of 'cbrtools'"
        for text in process.stderr:
            out_err.write(text)

        if process.wait() != 0:
            out_err.seek(0)
            error_str = out_err.read()
            raise Exception(error_str)

        out_string.seek(0)
        return DesignPrimersResults.from_json(json.load(out_string))

