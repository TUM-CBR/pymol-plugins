<<<<<<< HEAD
from typing import Callable, cast, Generic, NamedTuple, List, TypeVar, TYPE_CHECKING
=======
from typing import Callable, cast, NamedTuple, List, TypeVar
>>>>>>> main

class Point2d(NamedTuple):
    x : float
    y : float

TMeta = TypeVar('TMeta')

TNewMeta = TypeVar('TNewMeta')

# Python's type is a mess, therefore one can only have
# generics in namedtuples during typechecking.
if TYPE_CHECKING:
    class Series(NamedTuple, Generic[TMeta]):
        metadata : TMeta
        values : List[Point2d]

        def filter(self, fn: Callable[[Point2d], bool]) -> 'Series[TMeta]':
            return self._replace(
                values = list(filter(fn, self.values))
            )

        def update_meta(self, meta: TNewMeta) -> 'Series[TNewMeta]':
            return cast(
                'Series[TNewMeta]',
                self._replace(
                    metadata = meta
                )
            )
else:
    from typing import Any
    class Series(NamedTuple):
        metadata : Any
        values : List[Point2d]

        def filter(self, fn: Callable[[Point2d], bool]) -> 'Series[TMeta]':
            return self._replace(
                values = list(filter(fn, self.values))
            )

        def update_meta(self, meta: TNewMeta) -> 'Series[TNewMeta]':
            return cast(
                'Series[TNewMeta]',
                self._replace(
                    metadata = meta
                )
            )

class GlobalAttributes(NamedTuple):
    molar_extinction: float
    distance: float
    measurement_interval: float
    concentration_units: float

class RunMetadata(NamedTuple):
    concentration : float = 0.0

class KineticsRun(NamedTuple):
    run_metadata : RunMetadata
    data: List[float]

class KineticsRuns(NamedTuple):
    global_attributes : GlobalAttributes
    runs: List[KineticsRun]

class SubstrateInhibitionModel(NamedTuple):
    v_max : float = 0
    km : float = 0
    ksi : float = 0
    beta : float = 0

    def apply(self, s : float):
        num = (1 + self.beta*s/self.ksi)*s*self.v_max
        den = (1+s/self.ksi)*s + self.km

        return num / den

class EvalModelMetadata(NamedTuple):
    model_name: str

class FitModelResult(NamedTuple):
    model : SubstrateInhibitionModel