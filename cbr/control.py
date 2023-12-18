from typing import Dict, Iterable, Optional, Generic, Tuple, TypeVar

TVIter = TypeVar("TVIter")

def viter(*args : Optional[TVIter]) -> Iterable[TVIter]:
    arg : Optional[TVIter] = None
    for arg in args:
        if arg is not None:
            yield arg

TKey = TypeVar('TKey')
TValue = TypeVar('TValue')

class UpdateContextManager(Generic[TValue]):

    def __init__(
        self,
        collection : Dict[TKey, TValue],
        key : TKey,
        default : TValue
    ):

        self.__collection = collection
        self.__key = key
        self.__default = default
        self.__value = collection[key] if key in collection else default

    @property
    def value(self):
        return self.__value

    @value.setter
    def value(self, value: TValue):
        self.__value = value

    def __enter__(self, *args, **kwargs):
        return self

    def __exit__(self, t_exn, v_exn: Optional[Exception], tb_exn):

        if v_exn is None:
            self.__collection[self.__key] = self.value

update_key = UpdateContextManager

T = TypeVar('T')

#itertools.pairwise not in python 3.7
def pairwise(iter: Iterable[T]) -> Iterable[Tuple[T,T]]:
    last = None

    for item in iter:
        if last is not None:
            yield (last, item)

        last = item