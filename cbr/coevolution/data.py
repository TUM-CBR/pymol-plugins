from typing import Any, Dict, List, NamedTuple, Optional

from ..core import namedtuple
from ..extra.CbrExtraInteractiveHandler import MessageParser, MessageSerializer

class Scoring(NamedTuple):
    """
    Wether to residues are co-evolving or not is decided based
    on a score which constitutes many criteria. It is possible
    to weight each criteria differently by changing the paramters
    provided here.
    """

    occurrence_weight: float = 1
    exclusivity_weight: float = 1
    symmetry_weight: float = 1
    confidence_weight: float = 1
    confidence_treshold: float = 0.025
    confidence_center: float = 0.5
    confidence_concaveness: float = 8

class Query(NamedTuple):
    positions: List[int]
    max_results: int
    scoring: Scoring = Scoring()
    included_residues: Optional[Dict[int, List[str]]] = None

class InteractiveRequest(NamedTuple):
    query: Optional[Query]

def scoring_serialize(scoring: Scoring) -> Dict[Any, Any]:
    return {
        "symmetry_weight": scoring.symmetry_weight,
        "exclusivity_weight": scoring.exclusivity_weight,
        "confidence_weight": scoring.confidence_weight,
        "occurrence_weight": scoring.occurrence_weight,
        "confidence_treshold": scoring.confidence_treshold,
        "confidence_center": scoring.confidence_center,
        "confidence_concaveness": scoring.confidence_concaveness
    }

def query_serializer(query: Query) -> Dict[Any, Any]:
    return {
        "positions": query.positions,
        "max_results": query.max_results,
        "scoring": scoring_serialize(query.scoring),
        "included_residues": query.included_residues
    }

def interactive_request_serializer_implementation(value: InteractiveRequest) -> Dict[Any, Any]:

    result: Dict[Any, Any] = {}

    if value.query is not None:
        result["query"] = query_serializer(value.query)

    return result

interactive_request_serializer : MessageSerializer[InteractiveRequest] = interactive_request_serializer_implementation

class CoevolutionEntry(NamedTuple):
    residue_1: str
    residue_2: str
    score: float
    score_occurence: float
    score_exclusivity: float
    score_confidence: float
    score_symmetry: float

class CoevolutionPosition(NamedTuple):
    position: int
    by_position: Dict[int, CoevolutionEntry]

class CoevolutionResults(NamedTuple):
    positions: Dict[int, CoevolutionPosition]

class InteractiveResponse(NamedTuple):
    coevolution: Optional[CoevolutionResults]

def interactive_response_parser_implementation(value: Dict[Any, Any]) -> InteractiveResponse:
    print("starting the parser...")
    return namedtuple.parse(InteractiveResponse, value)

interactive_response_parser : MessageParser[InteractiveResponse] = interactive_response_parser_implementation