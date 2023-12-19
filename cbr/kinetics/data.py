from typing import NamedTuple, List

class GlobalAttributes(NamedTuple):
    molar_extinction: float
    distance: float
    measurement_interval: float

class RunMetadata(NamedTuple):
    concentration : float = 0.0

class KineticsRun(NamedTuple):
    run_metadata : RunMetadata
    data: List[float]

class KineticsRuns(NamedTuple):
    global_attributes : GlobalAttributes
    runs: List[KineticsRun]
