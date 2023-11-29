import json
from typing import Optional

from PyQt5.QtCore import QProcess

from ..extra.main import run_cbr_tools, run_cbr_tools_interactive
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

def create_cascade(
    db_file : str,
    spec_file : str,
    fasta_file : str,
    target_identity : float,
    email : str
) -> QProcess:

    return run_cbr_tools_interactive(
        [
            "cascades",
            "create",
            f"{target_identity}",
            email,
            db_file,
            fasta_file,
            spec_file
        ]
    )