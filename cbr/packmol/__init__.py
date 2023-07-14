from typing import Optional

from ..core.path import with_unique_prefix
from ..core.process import simple_execute

basic_input_template = \
    """
    tolerance 2.0
    filetype pdb
    output %s.pdb

    nloop 100

    structure %s.pdb
    number 1
    fixed 0. 0. 0 0. 0. 0.
    radius 1.5
    end structure
    """

def packmol_executable():

    #todo: winbugs
    return "packmol"

def pack_structure(
    input_structure : str,
    output_structure_filename : Optional[str] = None
    ):

    output_structure_filename = output_structure_filename or with_unique_prefix(input_structure, "packed")
    packmol_input = basic_input_template % (output_structure_filename, input_structure)

    with open(output_structure_filename, 'w') as output_structure:
        simple_execute(
            [packmol_executable()],
            packmol_input,
            output_structure
        )
