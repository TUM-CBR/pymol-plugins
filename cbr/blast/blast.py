from os import path
import os
import shutil

from ..core.Context import Context

class Blast:

    def __init__(self, _ctx: Context) -> None:

        own_folder = path.dirname(__file__)
        self.__blast_folder = path.join(
            own_folder,
            "resources",
            "bin"
        )

    def __get_command(self, cmd: str) -> str:

        if os.name == 'nt':
            return path.join(
                self.__blast_folder,
                f"{cmd}.exe"
            )
        
        cmd_path = shutil.which(cmd)

        if cmd_path is not None:
            return cmd_path
        
        raise Exception(
            f"The command '{cmd}' is not on your path."
        )
    
    @property
    def tblastn(self) -> str:
        return self.__get_command("tblastn")
    
    @property
    def makeblastdb(self) -> str:
        return self.__get_command("makeblastdb")