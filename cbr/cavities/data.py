from pymol import cmd
from typing import Any, Dict, List, NamedTuple, Optional

K_POINTS_ID = "points_id"
K_MIN_VOLUME = "min_volume"
K_MAX_VOLUME = "max_volume"

class FindCavitiesArgs(NamedTuple):
    points_id: str
    min_volume: int
    max_volume: int

    def to_json_dict(self) -> Dict[Any, Any]:
        return {
            K_POINTS_ID: self.points_id,
            K_MIN_VOLUME: self.min_volume,
            K_MAX_VOLUME: self.max_volume
        }
    

K_FIND_CAVITIES = 'find_cavities'

class CavitiesInteractiveInput(NamedTuple):
    """Class representing the input that is given to the cbr-extra-process
    which takes care of finding cavities in a protein."""
    find_cavities: Optional[List[FindCavitiesArgs]] = None

    def to_json_dict(self):

        v_find_cavities = self.find_cavities
        if v_find_cavities is not None:
            find_cavities = [
                find_cavity.to_json_dict()
                for find_cavity in v_find_cavities
            ]
        else:
            find_cavities = None

        return {
            K_FIND_CAVITIES: find_cavities
        }

    @classmethod
    def serialize(cls, value: 'CavitiesInteractiveInput') -> Dict[Any, Any]:
        return value.to_json_dict()

K_POINTS = 'points'
K_RADII = 'radii'

class CavityModel(NamedTuple):
    points: List[List[float]]
    radii: List[float]

    @classmethod
    def parse(cls, json_dict: Dict[Any, Any]) -> 'CavityModel':
        return CavityModel(
            points = json_dict[K_POINTS],
            radii = json_dict[K_RADII]
        )
    
    def volume(self) -> float:
        return sum(self.radii)
    
    def display(self, structure_name: str, chain: str, rep: str):
        for i,(pos,size) in enumerate(zip(self.points, self.radii)):
            cmd.pseudoatom(
                structure_name,
                pos=pos,
                elem='C',
                vdw=size/2,
                name='C1',
                resn='CAV',
                resi=str(i),
                chain=f"{chain}",
                segi='1'
            )
        
        cmd.show_as(
            representation=rep,
            selection = f"model {structure_name}"
        )

K_CAVITIES = "cavities"

class CavitiesResult(NamedTuple):
    cavities: Dict[str, List[CavityModel]]

    @classmethod
    def parse(cls, json_dict: Dict[Any, Any]) -> 'CavitiesResult':
        return CavitiesResult(
            cavities = {
                name : [CavityModel.parse(record) for record in records]
                for name, records in json_dict[K_CAVITIES].items()
            }
        )
    
    def get_cavity_name(self, name: str) -> str:
        name = name.replace("/", "_")
        return f"cav_{name}"
    
    def __display_cavity(self, name: str, rep: str):

        cavities = self.cavities[name]
        structure_name = self.get_cavity_name(name)
        cmd.delete(structure_name)

        for i,cavity in enumerate(cavities):
            cavity.display(structure_name, f"A{i}", rep)

    def display(self, rep: str):
        for cav in self.cavities:
            self.__display_cavity(cav, rep)

    
K_CAVITIES_RESULT = 'cavities_result'
    
class CavitiesInteractiveOutput(NamedTuple):
    """Class representing the output that is produced by the
    process which finds cavities."""
    cavities_result: Optional[CavitiesResult]

    @classmethod
    def parse(cls, json_dict: Dict[Any, Any]) -> 'CavitiesInteractiveOutput':

        cavities_json = json_dict.get(K_CAVITIES_RESULT)

        if cavities_json is not None:
            cavities_result = CavitiesResult.parse(cavities_json)
        else:
            cavities_result = None

        return CavitiesInteractiveOutput(cavities_result=cavities_result)
    
class AdvancedOptions(NamedTuple):
    empty_treshold: int = 0
    radii_scale: float = 1
    max_size_multiplier: float = 0.25

    def to_cmd_args(self):

        return [
            "--radii-scale",
            str(self.radii_scale),
            "--empty-treshold",
            str(self.empty_treshold),
            "--max_size_multiplier",
            str(self.max_size_multiplier)
        ]