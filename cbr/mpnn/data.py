from typing import cast, Dict, List, NamedTuple, Optional, Set, Tuple, Union
from pymol import cmd

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

    def add_residue(residue: str, resv: int):
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

        combined_exclusions : Dict[str, List[int]] = {}
        for position, res_list in self.excluded.items():
            residues = "".join(res_list)
            entry = combined_exclusions.get(residues)

            if entry is None:
                entry = combined_exclusions[residues] = []

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

class MpnnSpec(NamedTuple):
    edit_spaces: List[MpnnEditSpace]
    num_seqs : int
    mpnn_args: MpnnArgs

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
                residues = get_residues(model, chain)
                edit_positions = positions_to_be_edited[(model, chain)]
                fixed = [
                    residue.seq_pos + 1
                    for residue in residues.values()
                        if residue.resv not in edit_positions
                ]
                fixed.sort()
                results[model][chain] = fixed

        return results