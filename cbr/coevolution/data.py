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

    occurence_weight: float = 1
    exclusivity_weight: float = 1
    conserved_weight: float = 1

class Query(NamedTuple):
    positions: List[int]
    max_results: int
    scoring: Scoring = Scoring()

class InteractiveRequest(NamedTuple):
    query: Optional[Query]

def scoring_serialize(scoring: Scoring) -> Dict[Any, Any]:
    return {
        "occurence_weight": scoring.conserved_weight,
        "exclusivity_weight": scoring.exclusivity_weight,
        "conserved_weight": scoring.conserved_weight
    }

def query_serializer(query: Query) -> Dict[Any, Any]:
    return {
        "positions": query.positions,
        "max_results": query.max_results,
        "scoring": scoring_serialize(query.scoring)
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
    score_conserved: float

class CoevolutionPosition(NamedTuple):
    position: int
    by_position: Dict[int, CoevolutionEntry]

class CoevolutionResults(NamedTuple):
    positions: Dict[int, CoevolutionPosition]

class InteractiveResponse(NamedTuple):
    coevolution: Optional[CoevolutionResults]

def interactive_response_parser_implementation(value: Dict[Any, Any]) -> InteractiveResponse:
    return namedtuple.parse(InteractiveResponse, value)

interactive_response_parser : MessageParser[InteractiveResponse] = interactive_response_parser_implementation