from typing import cast, Dict, List, NamedTuple, Optional, Set, Tuple, Union
from pymol import cmd

from ..core.pymol.structure import StructureSelection

def mpnn_selection(model: str, chain : Optional[str] = None):
    model_selection = f"model {model}"

    if chain is not None:
        chain_selection = f" & chain {chain}"
    else:
        chain_selection = ""

    return model_selection + chain_selection + " & polymer.protein & backbone"

def get_chains(model: str) -> Set[str]:

    result: Set[str] = set()
    cmd.iterate(mpnn_selection(model), 'result.add(chain)', space={'result': result})
    return result

class Residue(NamedTuple):
    residue: str
    resv: int
    seq_pos: int

def get_residues(model: str, chain: str) -> Dict[int, Residue]:

    items : Dict[int, Residue] = {}
    valid_positions = get_valid_pdb_range(model, chain)

    def add_residue(residue: str, resv: int):

        if valid_positions.is_included(resv):
            items[resv] = Residue(residue=residue, resv=resv, seq_pos=0)

    cmd.iterate(
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


PositionsJonsl = Dict[str, Dict[str, List[int]]]
ExclusionListEntry = List[Union[List[int], str]]
ExclusionList = List[ExclusionListEntry]
ExclusionJonsl = Dict[str, Dict[str, ExclusionList]]

class MpnnEditSpace(NamedTuple):
    """
    Class representing a section of a protein to be generated using
    ProteinMPNN


    Attributes:
        model       The name of the model
        chain       The chain within the model
        residues    The residues that will be edited. This is the PDB value.
    """

    model: str
    chain: str
    residues: Set[int]
    excluded: Dict[int, Set[str]]

    def get_exclusion_jsonl(self) -> ExclusionJonsl:
        """
        Get a dict representing the jsonl that MPNN needs as an input for the exluded residues
        at particular positions. This jsonl looks like:
        {
            [model]: {
                [chain]: [
                    [position, "XXXXX"]
                ]
            }
        }
        """
        return {
            self.model : {
                self.chain : self.get_excluded_positions_array()
            }
        }

    def get_excluded_positions_array(self) -> ExclusionList:

        valid_range = get_valid_pdb_range(self.model, self.chain)
        combined_exclusions : Dict[str, List[int]] = {}
        for position, res_list in self.excluded.items():
            residues = "".join(res_list)
            entry = combined_exclusions.get(residues)

            if entry is None:
                entry = combined_exclusions[residues] = []

            if valid_range.is_included(position):
                entry.append(position)

        return [
            [positions, residues]
            for residues, positions in combined_exclusions.items()
        ]

class MpnnArgs(NamedTuple):
    """
    This class specifies the advanced options that can be passed to ProteinMPNN.

    Attributes:
        sampling_temperature    Sampling temperature for amino acids. Suggested values 0.1, 0.15, 0.2, 0.25, 0.3. Higher values will lead to more diversity.
        use_soluble_model       Flag to load ProteinMPNN weights trained on soluble proteins only.
        backbone_noise          Standard deviation of Gaussian noise to add to backbone atoms
    """

    sampling_temperature: float = 0.1
    use_soluble_model: bool = False
    backbone_noise: float = 0.0

TiedChainsEntryJonsl = Dict[str, List[int]]

TiedChainsJonsl = Dict[str, List[TiedChainsEntryJonsl]]

class TiedPositionsSpec(NamedTuple):
    tied_chains: List[Set[StructureSelection]]

    def __get_tied_positions(
        self,
        group: Set[StructureSelection],
        results_dict: TiedChainsJonsl
    ):
        items = list(group)

        if len(items) == 0:
            return
        
        model = items[0].structure_name
        assert all(i.structure_name == model and i.chain_name is not None for i in items), "Only chains in the same model can be tied"

        if model not in results_dict:
            results_dict[model] = []

        results_list = results_dict[model]
        chains = [cast(str, item.chain_name) for item in items]
        seqs = [get_residues(model, chain) for chain in chains]
        lengths = [len(s) for s in seqs]
        top = max(l for l in lengths)

        for i in range(top):
            tied_dict: TiedChainsEntryJonsl = {}
            for chain, length in zip(chains, lengths):
                if i >= length:
                    continue
                
                # The tied positions dict corresponds to all the values
                # that will be sampled together. That means that if we
                # want a position from multiple chains to be sampled together,
                # only that single position must be in the list per chain
                tied_dict[chain] = [i + 1]

            results_list.append(tied_dict)

        return results_list

    def get_tied_chains_jonsl(self) -> TiedChainsJonsl:

        results: TiedChainsJonsl = {}

        for group in self.tied_chains:
            self.__get_tied_positions(group, results)

        return results

class MpnnSpec(NamedTuple):
    edit_spaces: List[MpnnEditSpace]
    num_seqs : int
    mpnn_args: MpnnArgs
    tied_positions: TiedPositionsSpec

    def get_models(self) -> Set[str]:

        return set(
            space.model
            for space in self.edit_spaces
        )
    
    def get_excluded_jonsl(self) -> ExclusionJonsl:
        result : ExclusionJonsl = {}

        for space in self.edit_spaces:
            model = space.model
            chain = space.chain
            pdb_to_seq = get_residues(model, chain)

            chain_dict = result.get(model)
            if chain_dict is None:
                chain_dict = result[model] = {}

            def map_to_seq_position(value: ExclusionListEntry) -> ExclusionListEntry:
                [positions, seq] = value
                assert isinstance(positions, list), "First value is always an int"
                positions = [pdb_to_seq[pos].seq_pos + 1 for pos in positions]
                return [positions, seq]

            chain_dict[chain] = [
                map_to_seq_position(mapping)
                for mapping in space.get_excluded_positions_array()
            ]

        # ProteinMPNN requires that all models and chains have an entry
        # in the dictionary of excludsions. Hence we must add an empty
        # list for all models/chains that are not in any of the edit
        # spaces
        for model in self.get_models():
            if model not in result:
                result[model] = {}
            chains_dict = result[model]
            for chain in get_chains(model):
                if chain not in chains_dict:
                    chains_dict[chain] = []

        return result
    
    def get_chains_jsonl(self) -> Dict[str, List[List[str]]]:

        models = self.get_models()
        included : Dict[str, Set[str]] = dict((model, set()) for model in models)
        excluded : Dict[str, Set[str]] = dict((model, get_chains(model)) for model in models)

        for space in self.edit_spaces:
            included[space.model].add(space.chain)
            excluded[space.model].remove(space.chain)

        return {
            model: [list(included[model]), list(excluded[model])]
            for model in models
        }
    
    def get_positions_jsonl(self) -> PositionsJonsl:
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

        for model in self.get_models():
            results[model] = {}
            for chain in get_chains(model):
                position_range = get_valid_pdb_range(model, chain)
                residues = get_residues(model, chain)
                edit_positions: Optional[Set[int]] = positions_to_be_edited.get((model, chain))

                # If the model/chain combination is not present in the
                # 'edit_positions' dictionary, it means that no positions
                # are to be edited.
                if edit_positions is None:
                    edit_positions = set()
                fixed = [
                    residue.seq_pos + 1
                    for residue in residues.values()
                        if residue.resv not in edit_positions \
                            and position_range.is_included(residue.resv)
                ]
                fixed.sort()
                results[model][chain] = fixed

        return results
    
    def get_tied_jonsl(self) -> TiedChainsJonsl:
        return self.tied_positions.get_tied_chains_jonsl()