import pymol

from typing import NamedTuple, Optional

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

def get_selection_sequece(selection : str) -> str:
    result = []

    pymol.cmd.iterate(
        "%s & guide & alt +A" % selection,
        'result.append(oneletter)',
        space={'result': result}
    )
    return "".join(result)

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

def get_structure_offset(selection : StructureSelection) -> int:
    offset = [None]
    pymol.cmd.iterate(
        selection.selection,
        # if a protein has 999999999 residues, you probably
        # have bigger problems in life than your PhD thesis
        'offset[0] = min(resv, offset[0] or 999999999)',
        space={'offset' : offset, 'min' : min}
    )

    if offset[0] is None:
        raise Exception('The required structure "%s" has not been loaded' % selection.structure_name)

    return offset[0]