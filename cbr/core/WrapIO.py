from io import StringIO, TextIOBase
from typing import Any, Callable, Optional, cast, TextIO, TypeVar

NewStreamSpec = TypeVar('NewStreamSpec', str, Callable[[], TextIO])

def stream_builder(spec : NewStreamSpec) -> TextIO:

    if isinstance(spec, str):
        return open(spec, 'r')
    elif callable(spec):
        return spec()
    else:
        raise ValueError("Cannot build a stream from '%s'" % str(spec))

InitStreamSpec = TypeVar('InitStreamSpec', str, None, Callable[[TextIO], None])

def get_init_stream(spec : InitStreamSpec) -> 'Callable[[TextIO], None] | None':

    if(isinstance(spec, str)):
        data = spec

        def python_is_the_worse_language(stream : TextIO):
            stream.write(data)

        return python_is_the_worse_language

    else:
        return cast(Any, spec)

LegitStream = TypeVar('LegitStream', TextIO, TextIOBase)

class WrapIO(object):

    def __init__(
        self,
        stream : Optional[LegitStream] = None,
        open_stream : NewStreamSpec = StringIO,
        init : InitStreamSpec = None):

        self.__stream = cast(TextIO, stream or stream_builder(open_stream))
        self.__init = get_init_stream(init)
        self.__percolate = not stream

    @property
    def stream(self) -> TextIO:
        return self.__stream

    def __enter__(self, *args, **kwargs):

        if self.__percolate:
            self.__stream.__enter__(*args, **kwargs)

        if self.__init:
            self.__init(self.__stream)
            self.__stream.seek(0)

        return self

    def __exit__(self, *args, **kwargs):

        if self.__percolate:
            self.__stream.__exit__(*args, **kwargs)