from io import TextIOBase
import os
from typing import TypeVar

from .WrapIO import WrapIO

InputStream = TypeVar('InputStream', str, TextIOBase)

def get_input_stream(arg : InputStream) -> WrapIO:

    if isinstance(arg, TextIOBase):
        return WrapIO(stream=arg)

    elif isinstance(arg, str) and os.path.exists(arg):
        return WrapIO(open_stream=arg)
    else:
        return WrapIO(init=arg)