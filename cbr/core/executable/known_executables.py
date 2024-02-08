import os
from os import path
from PyQt5.QtWidgets import QDialog, QWidget

from .data import *
from .visual.KnownExecutable import KnownExecutableDialog

windows_known_executables = [

    KnownExecutableLocation(
        name = "ProteinMPNN",
        location = "https://github.com/TUM-CBR/cbr-tools-extra/releases/download/nightly-20240206/cbr-tools-extra-20240206-win64.zip",
        relative_path = "cbrtools/cbrtools.exe",
        known_executable = KnownExecutables.ProteinMPNN
    )
]

def check_executable(
    executable: Executable
):
    return path.exists(executable.location)

def find_windows_executable(
    widget: QWidget,
    executable: KnownExecutables,
    exe_dir: str
) -> Executable:
    
    resource = next(
        location
            for location in windows_known_executables
                if location.known_executable == executable
    )

    result = Executable(
        location = path.join(exe_dir, resource.relative_path),
        executable_type = ExecutableType.Binary
    )

    # Basic check for the time being
    if check_executable(result):
        return result
    
    dialog = KnownExecutableDialog(
        resource,
        exe_dir
    )
    dialog.show()

    download_result = dialog.exec()

    if download_result == QDialog.DialogCode.Accepted:
        if check_executable(result):
            return result
        else:
            raise Exception("Could not download the component for unknown reasons.")
    
    error = dialog.error()
    if error is not None:
        raise error
    
    raise Exception("Could not download the component for unknown reasons.")


unix_binaries = {
    KnownExecutables.ProteinMPNN: 'protein_mpnn_run.sh'
}

def find_executable(
    widget: QWidget,
    executable: KnownExecutables,
    exe_dir: str
) -> Executable:
    
    if os.name == 'nt':
        return find_windows_executable(widget, executable, exe_dir)
    
    binary = unix_binaries[executable]

    import shutil
    if shutil.which(binary):
        return Executable(
            location=binary,
            executable_type=ExecutableType.Binary
        )
    else:
        raise Exception(f"The executable '{binary}' is not in the system's path.")

