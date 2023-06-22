from io import StringIO, TextIOBase
from typing import Any, Callable, cast, TypeVar

NewStreamSpec = TypeVar('NewStreamSpec', str, Callable[[], TextIOBase])

def stream_builder(spec : NewStreamSpec) -> TextIOBase:

    if isinstance(spec, str):
        return open(spec, 'r')
    elif callable(spec):
        return spec()
    else:
        raise ValueError("Cannot build a stream from '%s'" % str(spec))

InitStreamSpec = TypeVar('InitStreamSpec', str, None, Callable[[TextIOBase], None])

def get_init_stream(spec : InitStreamSpec) -> 'Callable[[TextIOBase], None] | None':

    if(isinstance(spec, str)):
        data = spec

        def python_is_the_worse_language(stream : TextIOBase):
            stream.write(data)

        return python_is_the_worse_language

    else:
        return cast(Any, spec)

class WrapIO(object):

    def __init__(
        self,
        stream : 'TextIOBase | None' = None,
        open_stream : NewStreamSpec = StringIO,
        init : InitStreamSpec = None):

        self.__stream = stream or stream_builder(open_stream)
        self.__init = get_init_stream(init)
        self.__percolate = not stream

    @property
    def stream(self) -> TextIOBase:
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