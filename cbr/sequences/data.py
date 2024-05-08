from typing import NamedTuple, Optional, Sequence

class SequencesConfig(NamedTuple):
    selected_database: Optional[str] = None

    def to_json_dict(self):
        return {
            'selected_database': self.selected_database
        }

class QueryOrganism(NamedTuple):
    tax_id: int
    name: str

class QueryResult(NamedTuple):
    id: str
    query_id: str
    organism: Optional[QueryOrganism]
    identity: float
    query_sequence: str
    hit_sequence: str
    mid_line: str
    file_location: Optional[str]

class QueryResults(NamedTuple):
    results: Sequence[QueryResult]