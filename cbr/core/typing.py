from typing import cast, Optional, TypeVar

TResult = TypeVar('TResult')

def from_just(value : Optional[TResult]) -> TResult:
    return cast(TResult, value)