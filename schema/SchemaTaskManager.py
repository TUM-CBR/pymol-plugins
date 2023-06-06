import os
from PyQt5.QtCore import pyqtSignal, QObject
import re
import subprocess

from .raspp import schemacontacts
from .raspp import rasppcurve
from .SchemaResult import SchemaResult

seq_name_re = re.compile('>(?P<name>\w+)')

error_parse_sequence = \
    """Could not read the parent sequences. Please ensure the sequences follow the format:
    >SequenceName1
    [sequence]

    >SequenceName2
    [sequence]
    """

error_no_main_sequence = \
    """Could not find a sequence called 'Main'. This sequence will be aligned to the structure."""

class SchemaTaskManager(QObject):

    SCHEMA_RESULT_PREFIX = 'SCHEMA_results'
    __results_updated_signal = pyqtSignal(list)

    def __init__(self, schema_context, name : str, working_directory : str):
        super(SchemaTaskManager, self).__init__()
        self.__schema_context = schema_context
        self.__name = name
        self.__working_directory = working_directory

    def subscribe_results_updated(self, action):
        action(self.get_results())
        return self.__results_updated_signal.connect(action)

    def get_results(self):
        self.__ensure_directories()
        return \
            [SchemaResult(self.structure_name, self.pdb_file, os.path.join(self.location, file)) for file in os.listdir(self.location) \
            if file.startswith(SchemaTaskManager.SCHEMA_RESULT_PREFIX)]

    @property
    def location(self):
        return os.path.join(self.__working_directory, self.__name)

    @property
    def name(self):
        return self.__name
    
    @property
    def structure_name(self):
        return 'SCHEMA_%s' % self.name
    
    @property
    def pdb_file(self):
        self.__ensure_directories()
        return os.path.join(self.location, "%s.pdb" % self.structure_name)

    @property
    def pdb_aln_file(self):
        self.__ensure_directories()
        return os.path.join(self.location, 'SCHEMA_pdb_%s.aln' % self.name)

    @property
    def msa_aln_file(self):
        self.__ensure_directories()
        return os.path.join(self.location, 'SCHEMA_all_%s.aln' % self.name)

    @property
    def contacts_file(self):
        self.__ensure_directories()
        return os.path.join(self.location, 'SCHEMA_contacts_%s.txt' % self.name)

    def get_results_file_name(self, n : int):
        self.__ensure_directories()
        return os.path.join(
            self.location,
            "%s_%s_n%i" % (SchemaTaskManager.SCHEMA_RESULT_PREFIX, self.name, n)
            )
    
    def __ensure_directories(self):
        if not os.path.exists(self.location):
            os.makedirs(self.location)

    def save_resource(self, name : str, value : str):
        self.__ensure_directories()
        with open(os.path.join(self.location, name), 'w') as resource:
            resource.write(value)

    def load_resource(self, name : str):
        try:
            with open(os.path.join(self.location, name), 'r') as resource:
                return resource.read()
        except FileNotFoundError:
            return None

    @staticmethod
    def parse_sequences(sequences):

        results = {}
        current_key = None
        
        for line in filter(lambda x: x.strip() != "", sequences.splitlines()):

            if match := seq_name_re.match(line):
                current_key = match['name']
                results[current_key] = ""
            elif current_key:
                results[current_key] += line + "\n"
            else:
                raise ValueError(error_parse_sequence)

        return results

    def __run_clustal(self, outfile : str):
        return subprocess.Popen(
            [self.__schema_context.clustal, '-i', '-', '--force', '--outfmt=clustal', '-o', outfile],
            text=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
            
    def __align_parent(self, sequences : dict[str,str]):

        if 'Main' not in sequences:
            raise ValueError(error_no_main_sequence)

        if self.structure_name not in sequences:
            raise Exception('The sequence of structure "%s" is not in the list.' % self.structure_name)

        clustal = self.__run_clustal(self.pdb_aln_file)

        clustal.stdin.write('>Main\n%s\n>%s\n%s' % (sequences['Main'], self.structure_name, sequences[self.structure_name]))

        clustal.stdin.close()

        clustal.wait()


    def __align_sequences(self, sequences: dict[str,str]):

        clustal = self.__run_clustal(self.msa_aln_file)

        for (name, sequence) in sequences.items():

            # The sequence of the structure does not need to be included
            # in the multiple sequence alignment
            if name != self.structure_name:
                clustal.stdin.write('>%s\n%s\n\n' % (name, sequence))

        clustal.stdin.close()

        clustal.wait()

    def __run_schema_contacts(self):
        args = {
            schemacontacts.ARG_PDB_FILE: self.pdb_file,
            schemacontacts.ARG_PDB_ALIGNMENT_FILE: self.pdb_aln_file,
            schemacontacts.ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE: self.msa_aln_file,
            schemacontacts.ARG_OUTPUT_FILE: self.contacts_file
        }

        schemacontacts.main_impl(args)

    def __run_schema_curve(self, suffling_points: list[int], min_fragment_size: int):

        for n in suffling_points:
            
            args = {
                rasppcurve.ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE: self.msa_aln_file,
                rasppcurve.ARG_CONTACT_FILE: self.contacts_file,
                rasppcurve.ARG_NUM_CROSSOVERS: n,
                rasppcurve.ARG_OUTPUT_FILE: self.get_results_file_name(n),
                rasppcurve.ARG_MIN_FRAGMENT_SIZE: min_fragment_size
            }

            rasppcurve.main_impl(args)

    def run_schema(self, sequences_str: str, shuffling_points: list[int], min_fragment_size: int):
        sequences = SchemaTaskManager.parse_sequences(sequences_str)
        self.__align_parent(sequences)
        self.__align_sequences(sequences)
        self.__run_schema_contacts()
        self.__run_schema_curve([shuffling_points], min_fragment_size)
                
