from typing import Any, Callable, Tuple, cast, Generic, NamedTuple, List, TypeVar, TYPE_CHECKING

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
    measurement_interval: int
    concentration_units: float

class RunMetadata(NamedTuple):
    concentration : float = 0.001

class RunVelocityMetadata(NamedTuple):
    pass

class KineticsRun(NamedTuple):
    run_metadata : RunMetadata
    data: List[float]

class KineticsRuns(NamedTuple):
    global_attributes : GlobalAttributes
    runs: List[KineticsRun]

    def periods(self):
        return min(
                len(run.data)
                for run in self.runs
        )
    
    def concentrations(self):
        return [
            self.global_attributes.concentration_units*run.run_metadata.concentration
            for run in self.runs
        ]

class SubstrateInhibitionModel(NamedTuple):
    v_max : float = 0
    km : float = 0
    ksi : float = 0
    beta : float = 0

    def apply(self, s : float):
        num = (1 + self.beta*s/self.ksi)*s*self.v_max
        den = (1+s/self.ksi)*s + self.km

        return num / den
    
    def is_valid(self):
        for value in self:
            if cast(Any, value) is None:
                return False
            
        return True
    
    def scale(self, factor: float) -> 'SubstrateInhibitionModel':
        return SubstrateInhibitionModel(
            v_max = self.v_max*factor,
            km = self.km*factor,
            ksi = self.ksi*factor,

            # Beta is a ratio, so it does not make
            # sense to scale
            beta = self.beta
        )
    
SubstrateInhibitionModelFitRange = Tuple[SubstrateInhibitionModel, SubstrateInhibitionModel]

class EvalModelMetadata(NamedTuple):
    model_name: str

class SimulationSpec(NamedTuple):
    periods : int # Number of times the measurement is repeated
    interval : int # Seconds between measurement

    def times(self):
        return range(
            0,
            self.periods * self.interval,
            self.interval
        )

class SimulateModelMetadata(NamedTuple):
    model_name: str
    initial_concentration: float
    spec : SimulationSpec

class FitModelResult(NamedTuple):
    model : SubstrateInhibitionModel