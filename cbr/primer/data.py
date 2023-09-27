from Bio.Seq import Seq
import json
from typing import Any, Dict, List, NamedTuple

CODON_SIZE = 3

class PrimerResult(NamedTuple):
    left_primer: str
    tm_left: float
    right_primer: str
    tm_right: float
    inner_seq: str
    tm_all: float
    position: int
    amino_acid: str

    @property
    def c_left_primer(self):
        return str(Seq(self.left_primer).complement())

    @property
    def c_right_primer(self):
        return str(Seq(self.right_primer).complement())

    @property
    def c_inner(self):
        return str(Seq(self.inner_seq).complement())
    
    def to_json(self) -> str:
        return json.dumps(self._asdict())

    @classmethod
    def from_json(cls, json_str: str):
        return cls(**json.loads(json_str))

class DesignPrimersResults(NamedTuple):
    primers: List[PrimerResult]
    plasmid: str
    codon_mappings: Dict[str, str]
    
    def to_json(self) -> str:
        data = {
            'primers': [primer._asdict() for primer in self.primers],
            'plasmid': self.plasmid,
            'codon_mappings': self.codon_mappings
        }
        return json.dumps(data)

    @classmethod
    def from_json(cls, data: Dict[str, Any]):
        primers = [PrimerResult(**primer_data) for primer_data in data['primers']]
        return cls(primers=primers, plasmid=data['plasmid'], codon_mappings=data['codon_mappings'])

class DesignPrimersArgs(NamedTuple):
    sequence: str
    start: int
    codon_count: int
    min_length: int
    max_length: int
    organism: str
    
    def to_json(self) -> str:
        return json.dumps(self._asdict())

    @classmethod
    def from_json(cls, json_str: str):
        return cls(**json.loads(json_str))
