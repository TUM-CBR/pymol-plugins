import json
from typing import Optional

from ..extra.main import run_cbr_tools
from .data import *

def query_organisms(
    db_file : str,
    args : Optional[QueryCascadeArgs] = None
) -> QueryCascadeResult:

    if args is not None:
        extra = [
            str(args.max_identity_treshold),
        ] + \
        [
            f"{step.step_id},{step.policy}"
            for step in args.steps
        ]
    else:
        extra = []

    result = run_cbr_tools(
        [
            "cascades",
            "query",
            db_file
        ] + extra,
        out_handler=lambda out_stream: QueryCascadeResult.load(json.load(out_stream))
    )

    assert result is not None, "Bug in the code, result is unexpetedly None"
    return result