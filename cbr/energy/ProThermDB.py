import csv
import re
from typing import Iterable, List, NamedTuple, Optional
mutation_re = "(?P<original>[A-Za-z])(?P<position>\\d{1,3})(?P<new>[A-Za-z])"

class Mutation(NamedTuple):
    original_residue : str
    position : int
    new_residue : str

class ProThermEntry(NamedTuple):
    pdb_code : str
    mutation : Optional[Mutation]
    ddG : Optional[float]
    pH : str

def parse_mutation(mutation : str) -> Optional[Mutation]:
    item = re.match(mutation_re, mutation)

    if item:
        return Mutation(
            original_residue = item.group('original'),
            position = int(item.group('position')),
            new_residue = item.group('new')
        )

def read_ddg(item : str) -> Optional[float]:
    try:
        return float(item)
    except ValueError:
        return None

def read_entry(entry : dict) -> ProThermEntry:
    return ProThermEntry(
        pdb_code = entry['PDB_wild'],
        mutation = parse_mutation(entry['MUTATION']),
        ddG = read_ddg(entry["∆∆G_(kcal/mol)"]),
        pH = entry['pH']
    )

def read_entries(entries : Iterable[dict]) -> List[ProThermEntry]:
    return [read_entry(e) for e in entries]

def parse(location : str) -> List[ProThermEntry]:
    with open(location, 'r') as db:
        return read_entries(csv.DictReader(db, delimiter='\t'))
