from enum import Enum
from typing import List, NamedTuple

class CreateCascadeDatabaseArgs(NamedTuple):
    fasta_file : str

class QueryStepPolicy(Enum):
    keep = "keep"
    replace = "replace"
    any = "any"

class QueryCascadeStep(NamedTuple):
    step_id : int
    policy : QueryStepPolicy

class QueryCascadeArgs(NamedTuple):
    steps : List[QueryCascadeStep]
    max_identity_treshold : float

class OrganismResultEntry(NamedTuple):
    tax_id : int
    name : str

    @staticmethod
    def load(org_data: dict) -> 'OrganismResultEntry':
        return OrganismResultEntry(tax_id=org_data["tax_id"], name=org_data["name"])

class QueryCascadeResultStepEntry(NamedTuple):
    step_id : int
    step_name : str
    identity : float

    @staticmethod
    def load(step_data: dict) -> 'QueryCascadeResultStepEntry':
        return QueryCascadeResultStepEntry(
            step_id=step_data["step_id"],
            step_name=step_data["step_name"],
            identity=step_data["identity"]
        )

class QueryCascadeResultOrganismEntry(NamedTuple):
    organism : OrganismResultEntry
    steps : List[QueryCascadeResultStepEntry]

    @staticmethod
    def load(org: dict) -> 'QueryCascadeResultOrganismEntry':
        organism_entry = OrganismResultEntry.load(org['organism'])
        steps = [QueryCascadeResultStepEntry.load(step) for step in org["steps"]]
        return QueryCascadeResultOrganismEntry(organism=organism_entry, steps=steps)

class QueryCascadeResult(NamedTuple):
    organisms : List[QueryCascadeResultOrganismEntry]

    @staticmethod
    def load(input_dict: dict) -> 'QueryCascadeResult':
        organisms = [QueryCascadeResultOrganismEntry.load(org) for org in input_dict["organisms"]]
        return QueryCascadeResult(organisms=organisms)


