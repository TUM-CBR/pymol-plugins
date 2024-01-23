from typing import Iterable, Optional, TypeVar

TValue = TypeVar('TValue')

def assert_values(items: Iterable[Optional[TValue]]) -> Iterable[TValue]:
    for item in items:
        assert item is not None, "Item unexpectedly None"
        yield item