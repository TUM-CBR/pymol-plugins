import os

NOT_FOUND_ERROR="""Only windows binaries are currently provided with 'cbr-tools'. The feature you
are trying to use requires 'cbr-tools-extra'. Plese obtain a copy at https://github.com/TUM-CBR/cbr-tools-extra
and ensure it is in your PATH.
"""

def cbrtools_bin():

    if os.name == 'nt':
        return os.path.join(
            os.path.dirname(__file__),
            "resources",
            "cbrtools",
            "cbrtools.exe"
        )

    else:
        import shutil
        if shutil.which("cbrtools"):
            return "cbrtools"
        else:
            raise Exception(NOT_FOUND_ERROR)