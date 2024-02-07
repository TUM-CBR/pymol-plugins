from typing import cast, Dict, List, NamedTuple, Optional, Set, Tuple
from pymol import cmd

def mpnn_selection(model: str, chain : Optional[str] = None):
    model_selection = f"model {model}"

    if chain is not None:
        chain_selection = f" & chain {chain}"
    else:
        chain_selection = ""

    return model_selection + chain_selection

def get_chains(model: str) -> Set[str]:

    result: Set[str] = set()
    cmd.iterate(mpnn_selection(model), 'result.add(chain)')
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

class MpnnSpec(NamedTuple):
    edit_spaces: List[MpnnEditSpace]
    num_seqs : int

    def get_models(self) -> Set[str]:

        return set(
            space.model
            for space in self.edit_spaces
        )
    
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
                    residue.seq_pos
                    for residue in residues.values()
                        if residue.resv not in edit_positions
                ]
                fixed.sort()
                results[model][chain] = fixed

        return results