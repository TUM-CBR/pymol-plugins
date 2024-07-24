from typing import Any, Dict, List, NamedTuple, Optional, Sequence

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

class ErrorResult(NamedTuple):
    message: str
    offending_file: str

class ErrorResults(NamedTuple):
    results: Sequence[ErrorResult]

class SearchArg(NamedTuple):
    search_id: str
    tax_ids: Optional[List[int]] = None
    names: Optional[List[str]] = None
    accession: Optional[str] = None

    def to_json_dict(self) -> Dict[Any, Any]:
        return {
            'search_id': self.search_id,
            'tax_ids': self.tax_ids,
            'names': self.names,
            'accession': self.accession
        }

class SearchArgs(NamedTuple):
    searches: List[SearchArg]

    def to_json_dict(self) -> Dict[Any, Any]:
        return {
            'searches': [search.to_json_dict() for search in self.searches]
        }

class Organism(NamedTuple):
    taxid: int
    name: str

    @classmethod
    def from_json(cls, json_dict: Dict[str, Any]) -> 'Organism':
        return Organism(
            taxid = json_dict['taxid'],
            name = json_dict['name']
        )

    @classmethod
    def to_json(cls, organism: 'Organism') -> Dict[str, Any]:
        return {
            'taxid': organism.taxid,
            'name': organism.name
        }

class SearchResultRecord(NamedTuple):
    accession: str
    organism: Organism

    @classmethod
    def from_json(cls, json_dict: Dict[str, Any]) -> 'SearchResultRecord':
        return SearchResultRecord(
            accession = json_dict['accession'],
            organism = Organism.from_json(json_dict['organism'])
        )

    @classmethod
    def to_json(cls, record: 'SearchResultRecord') -> Dict[str, Any]:
        return {
            'accession': record.accession,
            'organism': Organism.to_json(record.organism)
        }

class SearchResultError(NamedTuple):
    message: str

    @classmethod
    def from_json(cls, json_dict: Dict[str, Any]) -> 'SearchResultError':
        return SearchResultError(
            message = json_dict['message']
        )

class SearchResult(NamedTuple):
    search_id: str
    records: List[SearchResultRecord]
    errors: List[SearchResultError]

    @classmethod
    def from_json(cls, json_dict: Dict[str, Any]) -> 'SearchResult':
        return SearchResult(
            search_id = json_dict['search_id'],
            records = [
                SearchResultRecord.from_json(record)
                for record in (json_dict['records'])
            ] if 'records' in json_dict else [],
            errors = [
                SearchResultError.from_json(error)
                for error in json_dict['errors']
            ] if 'errors' in json_dict else []
        )

class SaveSearchResult(NamedTuple):
    errors : Optional[List[str]] = None

    @classmethod
    def from_json(cls, json_dict: Dict[str, Any]) -> 'SaveSearchResult':
        return SaveSearchResult(
            errors = json_dict.get('errors')
        )

class InteractiveOutput(NamedTuple):
    search_result: Optional[SearchResult] = None
    save_search_result: Optional[SaveSearchResult] = None

    @classmethod
    def from_json_dict(cls, json_dict: Dict[str, Any]) -> 'InteractiveOutput':

        search_result = json_dict.get('search_result')
        save_search_result = json_dict.get('save_search_result')
        return InteractiveOutput(
            search_result = SearchResult.from_json(search_result) if search_result is not None else None,
            save_search_result = SaveSearchResult.from_json(save_search_result) if save_search_result is not None else None
        )

class InteractiveInput(NamedTuple):
    """Class representing the input that is provided to the
    "sequences" mode of the cbr-extra-process.
    """

    search: Optional[SearchArgs] = None
    save_search: Optional[List[SearchResultRecord]] = None

    def to_json_dict(self) -> Dict[Any, Any]:
        return {
            'search': self.search.to_json_dict() if self.search is not None else None,
            'save_search': [SearchResultRecord.to_json(record) for record in self.save_search] if self.save_search is not None else None
        }
