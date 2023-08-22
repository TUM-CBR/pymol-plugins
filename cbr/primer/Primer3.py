import os
import subprocess

class Primer3:

    @property
    def __oligotm_exe(self) -> str:

        if os.name == 'nt':
            os.path.join(
                os.path.dirname(__file__),
                "resources",
                "oligotm.exe"
            )

        return "oligotm"

    def __exception(self, text : str) -> Exception:
        return Exception(text)

    def oligo_tm(self, oligo : str) -> float:

        process = subprocess.Popen(
            [
                self.__oligotm_exe,
                oligo
            ],
            text = True,
            stdin = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )

        with process:

            if process.wait() != 0:
                raise self.__exception("Failed to launch 'oligotm' at location '%s'" % self.__oligotm_exe)

            if process.stdout:
                return float(process.stdout.read())

            raise self.__exception("Failed to open stdout for 'oligotm'")