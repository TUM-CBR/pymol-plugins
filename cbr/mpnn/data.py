from types import MappingProxyType
from typing import *
from pymol import cmd

from ..core.pymol.structure import StructureSelection

def mpnn_selection(
    model: str,
    chain : Optional[str] = None,
    include_ligands: bool = False
):
    model_selection = f"model {model}"

    if chain is not None:
        chain_selection = f" & chain {chain}"
    else:
        chain_selection = ""

    constraints = ["polymer.protein & backbone"]

    if include_ligands:
        constraints.append("organic")

    constraints_str = " | ".join(f"({constraint})" for constraint in constraints)

    return model_selection + chain_selection + f" & ({constraints_str})"

def get_chains(model: str) -> Set[str]:

    result: Set[str] = set()
    cmd.iterate(mpnn_selection(model), 'result.add(chain)', space={'result': result})
    return result

class Residue(NamedTuple):
    residue: str
    resv: int
    seq_pos: int

def get_residues(
    model: str,
    chain: str,
    state: Optional[int] = None
) -> Dict[int, Residue]:

    items : Dict[int, Residue] = {}
    valid_positions = get_valid_pdb_range(model, chain)

    def add_residue(residue: str, resv: int):

        if valid_positions.is_included(resv):
            items[resv] = Residue(residue=residue, resv=resv, seq_pos=0)

    if state is None:
        cmd.iterate(
            mpnn_selection(model, chain),
            'add_residue(oneletter, resv)',
            space={'add_residue': add_residue}
        )
    else:
        cmd.iterate_state(
            state,
            mpnn_selection(model, chain),
            'add_residue(oneletter, resv)',
            space={'add_residue': add_residue}
        )

    result = list(items.values())
    result.sort(key = lambda res: res.resv)
    return {
        res.resv: res._replace(seq_pos=seq_pos)
        for seq_pos, res in enumerate(result)
    }

class PdbRange(NamedTuple):
    min_value : int
    max_value : int

    def is_included(self, value: int) -> bool:
        return value >= self.min_value and value <= self.max_value

def get_valid_pdb_range(model: str, chain: str) -> PdbRange:
    """
    Gets the range of atoms that should be considered when providing
    positions to ProteinMPNN. The thing is that not all atoms end up
    in the pdb file, but ProteinMPNN will also add the missing positions
    if the missing atoms are in the middle.
    """
    selection = mpnn_selection(model, chain)
    min_value : Optional[int] = None
    max_value : Optional[int] = None

    def update_range(resv: int):
        nonlocal min_value
        nonlocal max_value
        if min_value is None or min_value > resv:
            min_value = resv
        if max_value is None or max_value < resv:
            max_value = resv

    cmd.iterate_state(0, selection, 'update_range(resv)', space={'update_range': update_range})

    if min_value is None or max_value is None:
        raise ValueError(f"There are no visible atoms in the model '{model}' and chain '{chain}'")

    return PdbRange(
        min_value=cast(int, min_value),
        max_value=cast(int, max_value)
    )


PositionsJonsl = Dict[str, List[str]]
ExclusionEntry = Dict[str, str]
ExclusionJonsl = Dict[str, ExclusionEntry]

class MpnnBackbone(NamedTuple):
    model : str
    chain : str

ResvToResidueMapping = MappingProxyType[int, Residue]

class MpnnEditSpace(NamedTuple):
    """
    Class representing a section of a protein to be generated using
    ProteinMPNN


    Attributes:
        backbone    The backbone correspnding to this edit space
        residues    The residues that will be edited. This is the PDB value.
    """
    backbone: MpnnBackbone
    residues: Set[int]
    excluded: Dict[int, Set[str]]

    @property
    def model(self) -> str:
        return self.backbone.model
    
    @property
    def chain(self) -> str:
        return self.backbone.chain

    def get_residues(self) -> MappingProxyType[int, Residue]:
        return MappingProxyType(get_residues(self.model, self.chain, state=1))

    def get_excluded_positions_array(self, pos_mapping: Callable[[int], int]) -> Iterable[Tuple[str, str]]:

        valid_range = get_valid_pdb_range(self.model, self.chain)
        combined_exclusions : Dict[str, List[int]] = {}
        for position, res_list in self.excluded.items():
            residues = "".join(res_list)
            entry = combined_exclusions.get(residues)

            if entry is None:
                entry = combined_exclusions[residues] = []

            if valid_range.is_included(position):
                entry.append(position)

        for residues, positions in combined_exclusions.items():
            for position in positions:
                yield (f"{self.chain}{pos_mapping(position)}", residues)

class MpnnArgs(NamedTuple):
    """
    This class specifies the advanced options that can be passed to ProteinMPNN.

    Attributes:
        sampling_temperature    Sampling temperature for amino acids. Suggested values 0.1, 0.15, 0.2, 0.25, 0.3. Higher values will lead to more diversity.
    """

    sampling_temperature: float = 0.1

TiedChainsJonsl = Dict[str, List[List[str]]]

class TiedPositionsSpec(NamedTuple):
    tied_chains: List[Set[StructureSelection]]

    def __get_tied_positions(
        self,
        group: Set[StructureSelection],
        results_dict: TiedChainsJonsl,
        models_to_pdb: Dict[str, str]
    ) -> None:
        items = list(group)

        if len(items) == 0:
            return
        
        model = items[0].structure_name
        pdb_path = models_to_pdb[model]
        assert all(i.structure_name == model and i.chain_name is not None for i in items), "Only chains in the same model can be tied"

        if pdb_path not in results_dict:
            results_dict[pdb_path] = []

        results_list = results_dict[pdb_path]
        chains_to_seqs = {
            chain: get_residues(model, chain, state=1)
            for item in items
            for chain in [cast(str, item.chain_name)]
        }

        for resv in set(resv for resvs in chains_to_seqs.values() for resv in resvs.keys()):

            symmetric_positions = [
                f"{chain}{resv}"
                for chain, seq in chains_to_seqs.items() if resv in seq
            ]

            if len(symmetric_positions) > 1:
                results_list.append(symmetric_positions)

    def get_tied_chains_jonsl(self, model_to_pdb: Dict[str, str]) -> TiedChainsJonsl:

        results: TiedChainsJonsl = {}

        for group in self.tied_chains:
            self.__get_tied_positions(group, results, model_to_pdb)

        return results

BackboneToResiduesMapping = MappingProxyType[MpnnBackbone, ResvToResidueMapping]

class MpnnSpec(NamedTuple):
    """Object which contains all the information needed to run one of the MPNN models."""

    edit_spaces: List[MpnnEditSpace]
    num_seqs : int
    mpnn_args: MpnnArgs
    tied_positions: TiedPositionsSpec
    mpnn_model: 'MpnnModel'

    def get_residues(self) -> BackboneToResiduesMapping:
        return MappingProxyType({
            space.backbone : space.get_residues()
            for space in self.edit_spaces
        })

    def get_models(self) -> Set[str]:

        return set(
            space.model
            for space in self.edit_spaces
        )
    
    def get_excluded_jonsl(self, model_to_pdb: Dict[str, str]) -> ExclusionJonsl:
        result : ExclusionJonsl = {}

        for space in self.edit_spaces:
            model_name = space.model
            model = model_to_pdb[model_name]

            entry = result.get(model)
            if entry is None:
                entry = result[model] = {}

            
            def position_mapping(pos: int) -> int:
                # return pdb_to_seq[pos].seq_pos + 1
                return pos
            
            for key,value in space.get_excluded_positions_array(position_mapping):
                entry[key] = value

        return result
    
    def get_chains_jsonl(self, model_to_pdb: Dict[str, str]) -> Dict[str, str]:

        models = self.get_models()
        included : Dict[str, Set[str]] = dict((model, set()) for model in models)
        excluded : Dict[str, Set[str]] = dict((model, get_chains(model)) for model in models)

        for space in self.edit_spaces:
            included[space.model].add(space.chain)
            excluded[space.model].remove(space.chain)

        return {
            model_to_pdb[model]: ",".join(included[model])
            for model in models
        }
    
    def get_positions_jsonl(self, model_to_pdb: Dict[str, str]) -> PositionsJonsl:
        """
        Get the dictionary of fixed positions for Protein MPNN. This dictionary is
        given relative to the residue position in the sequence, so we must convert
        the selected positions into sequence relative positions.
        """

        positions_to_be_edited : Dict[Tuple[str, str], Set[int]] = dict()

        for space in self.edit_spaces:
            key = (space.model, space.chain)
            positions = positions_to_be_edited.get(key)

            if positions is None:
                positions_to_be_edited[key] = positions = cast(Set[int], set())

            positions.update(space.residues)

        results : PositionsJonsl = {}

        for model_name in self.get_models():
            model = model_to_pdb[model_name]
            results[model] = []
            for chain in get_chains(model_name):
                position_range = get_valid_pdb_range(model_name, chain)
                residues = get_residues(model_name, chain)
                edit_positions: Optional[Set[int]] = positions_to_be_edited.get((model_name, chain))

                # If the model/chain combination is not present in the
                # 'edit_positions' dictionary, it means that no positions
                # are to be edited.
                if edit_positions is None:
                    edit_positions = set()
                fixed = [
                    f"{chain}{residue.resv}"
                    for residue in residues.values()
                        if residue.resv not in edit_positions \
                            and position_range.is_included(residue.resv)
                ]
                fixed.sort()
                results[model] += fixed

        return results
    
    def get_tied_jonsl(self, model_to_pdb: Dict[str, str]) -> TiedChainsJonsl:
        return self.tied_positions.get_tied_chains_jonsl(model_to_pdb)
    
class MpnnModel(NamedTuple):
    model_name: str
    include_ligand: bool = False

models : Dict[str, MpnnModel] = {
    "ProteinMPNN": MpnnModel("protein_mpnn"),
    "Soluble ProteinMPNN": MpnnModel("soluble_mpnn"),
    "LigandMPNN": MpnnModel("ligand_mpnn", include_ligand=True)
}
