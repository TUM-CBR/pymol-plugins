from io import StringIO, TextIOBase
from typing import Any, Callable, cast, TextIO, TypeVar

NewStreamSpec = TypeVar('NewStreamSpec', str, Callable[[], TextIO])

def get_stream_builder(spec : NewStreamSpec) -> Callable[[], TextIO]:

    if isinstance(spec, str):
        return lambda: open(spec, 'r')
    else:
        return cast(Callable, spec)

InitStreamSpec = TypeVar('InitStreamSpec', str, None, Callable[[TextIO], None])

def get_init_stream(spec : InitStreamSpec) -> 'Callable[[TextIO], None] | None':

    if(isinstance(spec, str)):
        data = spec

        def python_is_the_worse_language(stream):
            stream.write(data)

        return python_is_the_worse_language

    else:
        return cast(Any, spec)

class WrapIO(object):

    def __init__(
        self,
        stream : 'TextIO | None' = None,
        open_stream : NewStreamSpec = StringIO,
        init : InitStreamSpec = None):

        self.__stream = stream or get_stream_builder(open_stream).__call__()
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