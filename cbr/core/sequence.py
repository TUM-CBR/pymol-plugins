from Bio.Data.IUPACData import protein_letters_1to3

def residue_to_3(one : str) -> str:
    return protein_letters_1to3[one.upper()]