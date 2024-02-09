from Bio.Data.IUPACData import protein_letters_1to3, protein_letters_3to1, protein_letters

three_to_one_upper = dict(
    (three.upper(), one.upper())
    for three, one in protein_letters_3to1.items()
)

def residue_to_3(one : str) -> str:
    return protein_letters_1to3[one.upper()]

def residue_to_1(three: str) -> str:
    return three_to_one_upper[three.upper()]

def assert_residue(value: str) -> str:
    value = value.upper()
    assert value in protein_letters, f"The value '{value}' is not a known residue."
    return value