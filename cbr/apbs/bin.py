import os

def apbs_bin():

    if os.name == 'nt':
        return os.path.join(
            os.path.dirname(__file__),
            "resources",
            "APBS.windows",
            "bin",
            "apbs.exe"
        )

    return "apbs"

def pdb2pqr_bin():

    if os.name == 'nt':
        return os.path.join(
            os.path.dirname(__file__),
            "resources",
            "pdb2pqr.windows",
            "pdb2pqr.exe"
        )

    return "pdb2pqr"