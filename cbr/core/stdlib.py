from typing import TypeVar

T = TypeVar('T')

def get_or_raise(value : 'Exception | T') -> T:
    if isinstance(value, Exception):
        raise value
    else:
        return value