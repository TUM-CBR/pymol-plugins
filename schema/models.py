class SchemaResult(object):
    def __init__(self, pdb : str, results : list):
        self.__pdb = pdb
        self.__results = results

    @property
    def pdb(self):
        """The location of the pdb file that corresponds to these results"""
        return self.__pdb

    @property
    def results(self):
        """The results from running schema"""
        return self.__results
