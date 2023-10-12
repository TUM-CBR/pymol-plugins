from typing import Iterable, Optional, TypeVar

TVIter = TypeVar("TVIter")

def viter(*args : Optional[TVIter]) -> Iterable[TVIter]:
    arg : Optional[TVIter] = None
    for arg in args:
        if arg is not None:
            yield arg