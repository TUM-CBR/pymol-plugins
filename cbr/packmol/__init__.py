from os import path
import tempfile
from typing import Optional

from ..core.path import with_unique_prefix
from ..core.process import simple_execute

basic_input_template = \
    """
    tolerance 2.0
    filetype pdb
    output %s

    nloop 100

    structure %s
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

    # Packmol really really needs the input to be provided file redirect on the
    # stdin. It does very wierd stuff that no sane program shoud do.
    with tempfile.TemporaryDirectory() as packmol_workdir \
        , open(output_structure_filename, 'w') as output_structure:

        packmol_input_file_name = path.join(packmol_workdir, "packmol_inp.txt")

        with open(packmol_input_file_name, 'w') as packmol_inp_file:
            packmol_inp_file.write(packmol_input)

        simple_execute(
            [packmol_executable()],
            packmol_input_file_name,
            output_structure
        )
