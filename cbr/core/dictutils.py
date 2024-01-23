from typing import Callable, Dict, Iterable, Optional, TypeVar

TKey = TypeVar('TKey')
TValue = TypeVar('TValue')
TOut = TypeVar('TOut')

def merge(
    how: Callable[[Iterable[Optional[TValue]]], TOut],
    *dicts: Dict[TKey, TValue]) -> Dict[TKey, TOut]:

    keys = set(key for d in dicts for key in d.keys())

    return dict(
        (
            key,
            how(
                d.get(key)
                for d in dicts
            )
        )
        for key in keys
    )