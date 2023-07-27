AMINO_ACID_MAPPINGS = \
    { 'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K' \
    , 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N' \
    , 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W' \
    , 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M' }

INV_AMINO_ACID_MAPPING = dict((v,k) for (k, v) in AMINO_ACID_MAPPINGS.items())

def aa_from_letter(letter : str) -> str:
    return INV_AMINO_ACID_MAPPING[letter.upper()]