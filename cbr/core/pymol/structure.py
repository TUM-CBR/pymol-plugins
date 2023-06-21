import pymol

def get_structure_query(structure_name : str, chain : 'str | None' = None) -> str:
    if chain:
        return "(model %s) & (chain %s)" % (structure_name, chain)
    else:
        return "model %s" % structure_name

def get_pdb_sequence(structure_name : str, chain : 'str | None' = None) -> str:
    result = []

    pymol.cmd.iterate(
        "%s & guide & alt +A" % get_structure_query(structure_name, chain),
        'result.append(oneletter)',
        space={'result': result}
    )
    return "".join(result)

def get_structure_offset(structure_name : str, chain : 'str | None' = None) -> int:
    offset = [None]
    pymol.cmd.iterate(
        get_structure_query(structure_name, chain),
        # if a protein has 999999999 residues, you probably
        # have bigger problems in life than your PhD thesis
        'offset[0] = min(resv, offset[0] or 999999999)',
        space={'offset' : offset, 'min' : min}
    )

    if offset[0] is None:
        raise Exception('The required structure "%s" has not been loaded' % name)

    return offset[0] + 1