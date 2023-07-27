import json
import os
from os import path
import re
from typing import Iterable, NamedTuple, List, Optional

from ..gromacs.parsers import EnergyLogEntry

class Mutation(NamedTuple):
    selection : str
    mutation : str

    def to_json_dict(self):
        return {
            'selection' : self.selection,
            'mutation' : self.mutation
        }

    @staticmethod
    def from_json_dict(json_dict : dict):
        return Mutation(**json_dict)

class MutationContextBase(NamedTuple):

    name : str
    directory : str
    mutations : List[Mutation]

    def to_json_dict(self):
        return {
            'name' : self.name,
            'directory' : self.directory,
            'mutations': [m.to_json_dict() for m in self.mutations]
        }

    @staticmethod
    def from_json_dict(json_dict : dict):
        json_dict = json_dict.copy()
        K_MUTATIONS = "mutations"
        json_dict[K_MUTATIONS] = [Mutation.from_json_dict(m) for m in json_dict[K_MUTATIONS]]
        return MutationContextBase(**json_dict)

F_METADATA = "metadata.json"

class GromacsRun(NamedTuple):
    metadata : MutationContextBase
    emin_log : Optional[List[EnergyLogEntry]]

def parse_run(location : str) -> Optional[GromacsRun]:
    files = os.listdir(location)
    f_metadata = None
    f_emin_log = None

    for item in files:

        if item == F_METADATA:
            f_metadata = path.join(location, item)
        elif re.match("emin\\.\\w+\\.log", item):
            f_emin_log = path.join(location, item)

    if not f_metadata:
        return None

    with open(f_metadata, 'r') as stream:
        metadata = MutationContextBase.from_json_dict(json.load(stream))

    emin_log = None

    if f_emin_log:
        with open(f_emin_log, 'r') as stream:
            emin_log = EnergyLogEntry.parse_entries(stream)

    return GromacsRun(
        metadata = metadata,
        emin_log = emin_log
    )

def parse_runs(location : str) -> Iterable[GromacsRun]:

    runs = os.listdir(location)

    for run in runs:
        run = path.join(location, run)
        if path.isdir(run):
            result = parse_run(run)
            if result:
                yield result