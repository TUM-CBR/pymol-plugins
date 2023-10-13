from Bio.Seq import Seq
import csv
import json
from typing import Any, Dict, List, NamedTuple, TextIO

from ..core.Qt.visual.EditRecordsDialog import argsinfo

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

    def save_to_csv(self, out_stream : TextIO):

        csvwriter = csv.writer(out_stream)

        headers = [
            "Primer Position",
            "Codon",
            "Codon Residue",
            "Left Primer",
            "Tm Left Primer",
            "Right Primer",
            "Tm Right Primer",
            "Tm Amplicon"
        ]

        csvwriter.writerow(headers)
        csvwriter.writerows(
            [
                primer.position,
                primer.inner_seq,
                primer.amino_acid,
                primer.left_primer,
                primer.tm_left,
                primer.right_primer,
                primer.tm_right,
                primer.tm_all
            ]
            for primer in self.primers
        )

@argsinfo(
    mv_conc = 'Monovalent cation conc. (mM)',
    dv_conc = 'Divalent cation conc. (mM)',
    dntp_conc = 'dNTP conc. (mM)',
    dna_conc = 'DNA conc. (nM)',
    dmso_conc = "Concentration of DMSO (%)",
    dmso_fact = "DMSO correction factor, default 0.6",
    formamide_conc = "Concentration of formamide (mol/l)",
    annealing_temp_c = "Actual annealing temperature of the PCR reaction in (C)",
    max_nn_length = "Maximum length for nearest-neighbor calcs"
)
class Primer3Args(NamedTuple):
    mv_conc: float
    dv_conc: float
    dntp_conc : float
    dna_conc : float
    annealing_temp_c : float
    max_nn_length : int
    dmso_conc: float
    dmso_fact: float
    formamide_conc: float

    def to_json_dict(self) -> dict:

        json_dict = {}
        for field, field_type in self.__annotations__.items():
            assert field_type in [int, float], "Bug in code: This function needs to be updated!"
            json_dict[field] = getattr(self, field)

        return json_dict


DEFAULT_PRIMER3_ARGS = Primer3Args(
    mv_conc = 50.0,
    dv_conc = 1.5,
    dntp_conc = 0.2,
    dna_conc = 500,
    annealing_temp_c = -10.0,
    max_nn_length = 60,
    dmso_conc = 0.0,
    dmso_fact = 0.6,
    formamide_conc = 0.8
)

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
