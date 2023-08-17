import pymol

from typing import Dict, NamedTuple, Optional

class StructureSelection(NamedTuple):
    structure_name : str
    chain_name : Optional[str]
    segment_identifier : Optional[str]

    @property
    def selection(self) -> str:
        selectors = \
            [ self.structure_name
            , self.chain_name and "chain %s" % self.chain_name
            , self.segment_identifier and "segi %s" % self.segment_identifier
            ]

        return " and ".join(s for s in selectors if s)

def get_structure_query(structure_name : str, chain : 'str | None' = None) -> str:
    if chain:
        return "(model %s) & (chain %s)" % (structure_name, chain)
    else:
        return "model %s" % structure_name

def get_selection_sequence_index(selection : str) -> Dict[int, str]:
    result = []

    pymol.cmd.iterate(
        "%s & guide & alt +A" % selection,
        'result.append((int(resi), oneletter))',
        space={'result': result, 'int': int}
    )
    return dict(result)

def get_pdb_sequence_index(selection : StructureSelection) -> Dict[int, str]:
    return get_selection_sequence_index(
        selection.selection
    )

def get_selection_sequece(selection : str) -> str:
    sequence = get_selection_sequence_index(selection)
    return "".join(
        sequence[k]
        for k in sorted(sequence.keys())
    )

def get_pdb_sequence(selection : StructureSelection) -> str:
    return "".join(
        name
        for buffer in [[]] if pymol.cmd.iterate(
            "(%s) & guide & alt +A & bonded" % selection.selection,
            "buffer.append(oneletter)",
            space={'buffer': buffer}
        ) >= 0
        for name in buffer
    )