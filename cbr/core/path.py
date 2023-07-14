import itertools
from os import path

def with_unique_prefix(file_path : str, extra : str = "") -> str:
    directory = path.dirname(file_path)
    base_file_name = path.basename(file_path)

    for i in itertools.count(start=1):
        new_file_name = "%s%i.%s" % (extra, i, base_file_name)
        new_file = path.join(directory, new_file_name)
        if not path.exists(new_file):
            return new_file_name

    # Execution cannot reach this part of the code because the
    # iterator above should be infinite.
    raise Exception("Bug in the code!")