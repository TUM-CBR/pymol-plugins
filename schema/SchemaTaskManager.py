import os
from PyQt5.QtCore import pyqtSignal, pyqtSlot, QObject, QTimer
import re
import subprocess
from threading import Thread
from typing import Dict, List, NamedTuple

from .raspp import schemacontacts
from .raspp import rasppcurve
from .SchemaResult import SchemaResult
from .support import fasta

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

    __no_clustal_stdin = Exception("Bug in the code! The clustal process has no stdin")
    SCHEMA_RESULT_PREFIX = 'SCHEMA_results'
    __results_updated_signal = pyqtSignal(list)
    is_busy_signal = pyqtSignal(int)

    def __init__(self, schema_context):
        super(SchemaTaskManager, self).__init__()
        self.__schema_context = schema_context
        self.__current_tasks = []
        self.__schema_watcher = QTimer()
        self.__schema_watcher.setInterval(100)
        self.__schema_watcher.setSingleShot(False)
        self.__schema_watcher.timeout.connect(self.__check_tasks)
        self.__schema_watcher.start()
        self.__current_task_count = 0

    @property
    def __working_directory(self) -> str:
        return self.__schema_context.working_directory

    def __check_tasks(self):
        prev_task_count = self.__current_task_count
        self.__current_tasks = list(filter(lambda x: not x.is_done, self.__current_tasks))
        self.__current_task_count = len(self.__current_tasks)

        if prev_task_count == self.__current_task_count:
            return
        else:
            self.is_busy_signal.emit(self.__current_task_count)
            self.__results_updated_signal.emit(self.get_results())

    def subscribe_results_updated(self, action):
        action(self.get_results())
        return self.__results_updated_signal.connect(action)

    def get_pdb_file_name(self, name : str) -> str:
        return SchemaTask.get_pdb_file_name(self.__working_directory, name)

    def __get_results(self):

        try:
            for folder_name in os.listdir(self.__working_directory):
                folder = os.path.join(self.__working_directory, folder_name)
                try:
                    for file_name in os.listdir(folder):
                        file = os.path.join(folder, file_name)
                        if  file_name.startswith(SchemaTaskManager.SCHEMA_RESULT_PREFIX):
                            yield SchemaResult(
                                folder_name,
                                SchemaTask.get_structure_name(folder_name),
                                SchemaTask.get_pdb_file_name(self.__working_directory, folder_name),
                                file)
                except NotADirectoryError:
                    pass
        except FileNotFoundError:
            pass

    def __ensure_directories(self):
        if not os.path.exists(self.__working_directory):
            os.makedirs(self.__working_directory)

    def save_resource(self, name : str, value : str):
        self.__ensure_directories()
        with open(os.path.join(self.__working_directory, name), 'w') as resource:
            resource.write(value)

    def load_resource(self, name : str):
        try:
            with open(os.path.join(self.__working_directory, name), 'r') as resource:
                return resource.read()
        except FileNotFoundError:
            return None

    def get_results(self):
        return list(self.__get_results())

    def run_schema(self, name: str, sequences_str: str, shuffling_points: List[int], min_fragment_size: int) -> None:
        task = SchemaTask(self.__schema_context, name, self.__working_directory, sequences_str, shuffling_points, min_fragment_size)
        self.__current_tasks.append(task)
        self.__check_tasks()

class SchemaTask(QObject):

    is_busy_signal = pyqtSignal(bool)

    def __init__(
        self,
        schema_context,
        name : str,
        working_directory : str,
        sequences_str : str,
        shuffling_points : List[int],
        min_fragment_size: int):

        super(SchemaTask, self).__init__()
        self.__schema_context = schema_context
        self.__name = name
        self.__working_directory = working_directory
        self.__schema_thread = Thread(
            target = self.__run_schema_action,
            args = [sequences_str, shuffling_points, min_fragment_size]
        )

        self.__schema_thread.start()

    @property
    def is_done(self):
        return not self.__schema_thread.is_alive()

    @staticmethod
    def __location(working_directory : str, name : str):
        return os.path.join(working_directory, name)

    @property
    def location(self):
        return SchemaTask.__location(self.__working_directory, self.__name)

    @property
    def name(self):
        return self.__name
    
    @property
    def structure_name(self):
        return SchemaTask.get_structure_name(self.name)
    
    @property
    def pdb_file(self):
        self.__ensure_directories()
        return SchemaTask.get_pdb_file_name(self.__working_directory, self.name)

    @staticmethod
    def get_structure_name(name: str):
        return 'SCHEMA_%s' % name

    @staticmethod
    def get_pdb_aln_file_name(working_directory : str, result_name: str):
        return os.path.join(
            working_directory,
            result_name,
            'SCHEMA_pdb_%s.aln' % result_name
        )

    @staticmethod
    def get_pdb_file_name(working_directory : str, result_name: str):
        SchemaTask.ensure_location(SchemaTask.__location(working_directory, result_name))
        return os.path.join(
            working_directory,
            result_name,
            "%s.pdb" % SchemaTask.get_structure_name(result_name)
        )

    @property
    def pdb_aln_file(self):
        self.__ensure_directories()
        return SchemaTask.get_pdb_aln_file_name(self.__working_directory, self.name)

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
    
    @staticmethod
    def ensure_location(location : str) -> None:
        if not os.path.exists(location):
            os.makedirs(location)

    def __ensure_directories(self):
        SchemaTask.ensure_location(self.location)

    @staticmethod
    def parse_sequences(sequences : str) -> Dict[str,str]:

        results : Dict[str,str] = {}
        
        for line in fasta.parse_fasta_iter(sequences):

            if isinstance(line, BaseException):
                raise line
            else:
                (k,v) = line
                results[k] = v

        return results

    def __run_clustal(self, outfile : str) -> subprocess.Popen:
        return subprocess.Popen(
            [self.__schema_context.clustal, '-i', '-', '--force', '--outfmt=clustal', '-o', outfile],
            text=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
            
    def __align_parent(self, sequences : Dict[str,str]) -> None:

        if 'Main' not in sequences:
            raise ValueError(error_no_main_sequence)

        if self.structure_name not in sequences:
            raise Exception('The sequence of structure "%s" is not in the list.' % self.structure_name)

        with self.__run_clustal(self.pdb_aln_file) as clustal:
            stdin = clustal.stdin
            if stdin:
                stdin.write('>Main\n%s\n>%s\n%s' % (sequences['Main'], self.structure_name, sequences[self.structure_name]))
                stdin.close()
            else:
                raise SchemaTaskManager.__no_clustal_stdin

            clustal.wait()


    def __align_sequences(self, sequences: Dict[str,str]):

        with self.__run_clustal(self.msa_aln_file) as clustal:
            if not clustal.stdin:
                raise SchemaTaskManager.__no_clustal_stdin

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

    def __run_schema_curve(self, suffling_points: List[int], min_fragment_size: int):

        for n in suffling_points:
            
            args = {
                rasppcurve.ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE: self.msa_aln_file,
                rasppcurve.ARG_CONTACT_FILE: self.contacts_file,
                rasppcurve.ARG_NUM_CROSSOVERS: n,
                rasppcurve.ARG_OUTPUT_FILE: self.get_results_file_name(n),
                rasppcurve.ARG_MIN_FRAGMENT_SIZE: min_fragment_size
            }

            rasppcurve.main_impl(args)

    def __run_schema_action(self, sequences_str: str, shuffling_points: List[int], min_fragment_size: int):
        sequences = SchemaTask.parse_sequences(sequences_str)
        self.__align_parent(sequences)
        self.__align_sequences(sequences)
        self.__run_schema_contacts()
        self.__run_schema_curve(shuffling_points, min_fragment_size)
