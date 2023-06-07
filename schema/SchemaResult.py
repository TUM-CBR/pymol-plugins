import os

class SchemaResultEntry(object):

    def __init__(self, line : str):
        entry = line.replace("\t", " ").split()
        self.__energy = entry[0]
        self.__mutations = entry[1]
        self.__shuffling_points = [int(x) for x in entry[2:]]

    @property
    def energy(self):
        return self.__energy

    @property
    def mutations(self):
        return self.__mutations

    @property
    def shuffling_points(self):
        return self.__shuffling_points

class SchemaResult(object):
    def __init__(
        self,
        result_name : str,
        structure_name : str,
        pdb : str,
        results : str):
        
        self.__result_name = result_name
        self.__pdb = pdb
        self.__results = results
        self.__structure_name = structure_name

    @property
    def structure_name(self):
        """Return the name of the structure corresponding to this result"""
        return self.__structure_name

    @property
    def pdb(self):
        """The location of the pdb file that corresponds to these results"""
        return self.__pdb

    @property
    def results_file(self):
        """The results from running schema"""
        return self.__results

    @property
    def name(self):
        return os.path.join(self.__result_name, os.path.basename(self.__results))

    def load_results(self):

        return_value = []
        with open(self.__results, 'r') as results:
            for result in results:

                # The line is a comment, igonre
                if result.startswith('#'):
                    continue
                else:
                    return_value.append(SchemaResultEntry(result))

        return return_value