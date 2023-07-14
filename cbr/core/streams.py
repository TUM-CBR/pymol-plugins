import io
import os
from typing import TextIO, TypeVar

from .WrapIO import WrapIO

InputStream = TypeVar('InputStream', str, TextIO, io.TextIOBase)

def get_input_stream(arg : InputStream) -> WrapIO:

    if isinstance(arg, io.TextIOBase):
        return WrapIO(stream=arg)
    elif isinstance(arg, str) and os.path.exists(arg):
        return WrapIO(open_stream=arg)
    elif isinstance(arg, str):
        return WrapIO(init=arg)
    else:
        raise ValueError("The 'arg' must be a string or a stream.")