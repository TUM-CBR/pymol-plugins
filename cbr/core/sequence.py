AMINO_ACID_MAPPINGS = \
    { 'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K' \
    , 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N' \
    , 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W' \
    , 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M' }

INV_AMINO_ACID_MAPPING = dict((v,k) for (k, v) in AMINO_ACID_MAPPINGS.items())

def aa_from_letter(letter : str) -> str:
    return INV_AMINO_ACID_MAPPING[letter.upper()]

def letter_from_aa(resi : str) -> str:
    return AMINO_ACID_MAPPINGS[resi.upper()]

def ensure_oneletter(resi : str) -> str:

    if resi in INV_AMINO_ACID_MAPPING:
        return resi
    else:
        return letter_from_aa(resi)

def ensure_abbreviation(resi : str) -> str:

    if resi in AMINO_ACID_MAPPINGS:
        return resi
    else:
        return aa_from_letter(resi)