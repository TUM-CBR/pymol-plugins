import csv
import glob
import os
from os import path
from typing import Iterable, List, NamedTuple, Optional, TextIO

from ..gromacs.parsers import EnergyLogEntry

from .mutations import MutationResult

class EnergyMinimizationRun(NamedTuple):
    result_metadata : MutationResult
    result_energy_log : List[EnergyLogEntry]

    @staticmethod
    def from_file(result_file : str) -> 'EnergyMinimizationRun':

        with open(result_file, 'r') as result_stream:
            result_metadata = MutationResult.from_stream(result_stream)

        energy_log_file = result_metadata.get_energy_log_file(path.dirname(result_file))
        with open(energy_log_file, 'r') as energy_log_stream:
            result_energy_log = EnergyLogEntry.parse_entries(energy_log_stream)

        return EnergyMinimizationRun(
            result_metadata = result_metadata,
            result_energy_log = result_energy_log
        )

    @staticmethod
    def from_folder(result_folder : str) -> 'Optional[EnergyMinimizationRun]':
        metadata_file = next(glob.glob(path.join(result_folder, "*.metadata.json")).__iter__(), None)
        if metadata_file:
            return EnergyMinimizationRun.from_file(metadata_file)

    @staticmethod
    def from_folders(folder : str) -> 'Iterable[EnergyMinimizationRun]':

        for filename in os.listdir(folder):

            full_file_path = path.join(folder, filename)
            result = EnergyMinimizationRun.from_folder(full_file_path)
            if result:
                yield result


    def to_csv_entry(self, stream: TextIO):
        pass