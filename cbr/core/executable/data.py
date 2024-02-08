from enum import Enum
from typing import NamedTuple

class ExecutableType(Enum):
    Binary = 0

class Executable(NamedTuple):
    executable_type : ExecutableType
    location : str

class KnownExecutables(Enum):
    ProteinMPNN = "ProteinMPNN"

class KnownExecutableLocation(NamedTuple):
    name : str
    location: str
    relative_path: str
    known_executable: KnownExecutables
