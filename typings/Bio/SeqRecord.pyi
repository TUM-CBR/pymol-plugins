from typing import Any, Dict, List, Optional, Sequence, Union

from Bio.SeqFeature import SeqFeature
from Bio.Seq import MutableSeq, Seq

_AnnotationsDictValue = Union[str, int]
_AnnotationsDict = Dict[str, _AnnotationsDictValue]

class SeqRecord:
    _id: str
    _seq: Seq

    def __init__(
        self,
        seq: Optional[Union[Seq, MutableSeq]],
        id: Optional[str] = "<unknown id>",
        name: str = "<unknown name>",
        description: str = "<unknown description>",
        dbxrefs: Optional[List[str]] = None,
        features: Optional[List["SeqFeature"]] = None,
        annotations: Optional[_AnnotationsDict] = None,
        letter_annotations: Optional[Dict[str, Sequence[Any]]] = None,
    ) -> None: ...