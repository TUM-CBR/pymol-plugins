from typing import Generic, NamedTuple, List, TypeVar

class Point2d(NamedTuple):
    x : float
    y : float

TMeta = TypeVar('TMeta')

class Series(NamedTuple, Generic[TMeta]):
    metadata : TMeta
    values : List[Point2d]

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
