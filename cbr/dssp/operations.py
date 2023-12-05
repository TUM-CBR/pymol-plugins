from io import StringIO
import os
import subprocess

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
        if shutil.which("mkdssp"):
            return "mkdssp"
        else:
            raise Exception(NOT_FOUND_ERROR)


def run_dssp(
    input_file: str,
    output_file: str
) -> StringIO:

    process = subprocess.Popen(
        [dssp_exe(), input_file, output_file],
        text = True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    with process, \
        StringIO() as out_stream, \
        StringIO() as error_stream:

        assert process.stdout, f"Error opening the stdout of {dssp_exe()}"
        assert process.stderr, f"Error opening the stderr of {dssp_exe()}"

        for data in process.stdout:
            out_stream.write(data)

        for data in process.stderr:
            error_stream.write(data)

        if process.wait() != 0:
            error_stream.seek(0)
            raise Exception(error_stream.read())

        out_stream.seek(0)
        return out_stream
