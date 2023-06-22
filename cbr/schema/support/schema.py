from typing import Dict

from ...clustal import msa
from .strings import MAIN_SEQUENCE_NAME

class SchemaAlignments(object):

    @staticmethod
    def from_files(pdb_msa_file_name : str, parents_msa_file_name : str) -> 'SchemaAlignments':

        with open(pdb_msa_file_name, 'r') as pdb_msa_file:
            pdb_msa = msa.parse_alignments(pdb_msa_file.read())

        with open(parents_msa_file_name, 'r') as parents_msa_file:
            parents_msa = msa.parse_alignments(parents_msa_file.read())

        return SchemaAlignments(
            pdb_msa,
            parents_msa
        )

    def __init__(self, structure_msa : Dict[str, str], parents_msa : Dict[str, str]):

        if len(structure_msa) != 2:
            raise ValueError("Only one parent must be aligned to the structure")

        if MAIN_SEQUENCE_NAME not in structure_msa \
            or MAIN_SEQUENCE_NAME not in parents_msa:
            raise ValueError("A sequence called %s must be in both msa" % MAIN_SEQUENCE_NAME)

        for (name, sequence) in structure_msa.items():

            if name == MAIN_SEQUENCE_NAME:
                self.__main_msa = sequence
            else:
                self.__structure_msa = sequence

        self.__parents_msa = parents_msa

    @property
    def __main_parents_msa(self) -> str:
        return self.__parents_msa[MAIN_SEQUENCE_NAME]

    def get_structure_position(self, i : int) -> int:
        """
        This function converts an index relative to a position in the msa file of all
        parents to the position that index corresponds in the structure.
        """

        main_parents_prefix = msa.remove_spacers(self.__main_parents_msa[0:i])
        main_parents_prefix_pos = 0
        result = 0

        while(main_parents_prefix_pos < len(main_parents_prefix)):
            parent_current = main_parents_prefix[main_parents_prefix_pos]
            seq_current = self.__main_msa[result]
            result += 1

            if seq_current == parent_current:
                main_parents_prefix_pos += 1
            elif not msa.is_spacer(seq_current):
                raise ValueError("The character '%s' is not a valid clustal character." % seq_current)

        return len(msa.remove_spacers(self.__structure_msa[0:result - 1])) - 1

