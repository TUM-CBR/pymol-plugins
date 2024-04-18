from typing import Any, Dict, List, NamedTuple, Optional, Type, TypeVar

from ..extra.CbrExtraInteractiveHandler import CbrExtraInteractiveHandler, MessageParser, MessageSerializer

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

def coevolution_entry_parser(value: Any) -> CoevolutionEntry:

    if not isinstance(value, dict):
        raise ValueError("value: Expected a dictionary")
    
    return CoevolutionEntry(
        residue_1=value["residue_1"],
        residue_2=value["residue_2"],
        score=value["score"],
        score_occurence=value["score_occurence"],
        score_exclusivity=value["score_exclusivity"],
        score_conserved=value["score_exclusivity"]
    )

def coevolution_position_parser(value: Any) -> CoevolutionPosition:

    if not isinstance(value, dict):
        raise ValueError("value: Expected a dictionary")
    
    position = value.get("position")
    by_position = value.get("by_position")

    if not isinstance(position, int):
        raise ValueError("value: Expected a position key with an int value")
    
    if not isinstance(by_position, dict):
        raise ValueError("value: Expected a by_position key with a dict value")
    
    return CoevolutionPosition(
        position=position,
        by_position={
            k: coevolution_entry_parser(v)
            for k,v in by_position.items()
        }
    )

def coevolution_parser(value: Any) -> CoevolutionResults:

    if not isinstance(value, dict):
        raise ValueError("value: Expected a dictionary")
    
    positions: Any = value.get("positions")

    if not isinstance(positions, dict):
        raise ValueError("value: Expected to have a 'positions' key.")

    return CoevolutionResults(
        positions={
            k: coevolution_position_parser(v)
            for k,v in positions.items()
        }
    )


def interactive_response_parser_implementation(value: Dict[Any, Any]) -> InteractiveResponse:

    coevolution = value.get('coevolution')

    if coevolution is not None:
        return InteractiveResponse(
            coevolution=
        )
    else:
        raise ValueError("value: Unknown message")

interactive_response_parser : MessageParser[InteractiveResponse] = interactive_response_parser_implementation