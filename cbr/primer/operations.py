from io import StringIO
import json

from ..extra.main import run_cbr_tools
from .data import DesignPrimersResults, Primer3Args

def design_primers(
    start : int,
    codon_count : int,
    min_length : int,
    max_length : int,
    organism : str,
    sequence : str,
    primer3Args : Primer3Args,
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
                "sequence": sequence,
                "primer3Args": primer3Args.to_json_dict()
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