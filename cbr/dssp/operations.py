import os

NOT_FOUND_ERROR="""The executable "mkdssp" is needed to use this feature."""

def dssp_exe():

    if os.name == 'nt':
        return os.path.join(
            os.path.dirname(__file__),
            "resources",
            "mkdssp.exe"
        )

    else:
        import shutil
        if shutil.which("cbrtools"):
            return "cbrtools"
        else:
            raise Exception(NOT_FOUND_ERROR)