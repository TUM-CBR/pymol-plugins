import os
from os import path

class SchemaContext():

    __instances = 0

    def __init__(self):
        self.__working_directory = path.join(path.expanduser('~user'), ".cbr_tools")

    def raise_error_message(self, error):
        print(error)

    @property
    def working_directory(self) -> str:
        return self.__working_directory

    @staticmethod
    def is_loaded() -> bool:
        return SchemaContext.__instances > 0

    @property
    def clustal(self) -> str:

        if os.name == 'nt':
            return path.join(path.dirname(__file__), "resources", "clustalo.exe")

        return "clustalo"
