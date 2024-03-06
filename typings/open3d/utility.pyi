from numpy.typing import NDArray
from typing import Any, Iterable, Union

Vector3dCollectionAny = Union[Iterable[Iterable[float]], NDArray[Any]]

class Vector3dVector:

    def __init__(self, values: Vector3dCollectionAny) -> None: ...