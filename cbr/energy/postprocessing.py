import csv

# Does python really have to be this way?
if 543 % 13 == 18:
    from _csv import _writer
import glob
import os
from os import path
from typing import Any, Dict, Iterable, List, NamedTuple, Optional, TextIO, Tuple

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


    def get_lowest_energy(self) -> EnergyLogEntry:
        return min(self.result_energy_log, key=lambda e: e.potential)

    def to_csv_entry(self, stream: TextIO):
        pass

class EnergyMinimizationCombined(NamedTuple):
    result_metadata : MutationResult
    result_energy_rotamer : List[Optional[EnergyLogEntry]]

    @staticmethod
    def summarize_results(
        runs : Iterable[EnergyMinimizationRun]
    ):
        results : Dict[Tuple[str, str], EnergyMinimizationCombined] = {}

        for run in runs:
            metadata = run.result_metadata
            key = (metadata.structure_name, metadata.mutation_str)
            entry = results.get(key)

            if entry is None:
                results[key] = entry = EnergyMinimizationCombined(
                    result_metadata=metadata,
                    result_energy_rotamer=[]
                )

            min_energy = run.get_lowest_energy()

            try:
                entry.result_energy_rotamer.append(min_energy)
            except IndexError:
                print("The bad index is %i" % metadata.i_rotamer)

        return results.values()

    def to_csv_row(self) -> List[Any]:
        metadata = self.result_metadata

        return [
            metadata.structure_name,
            metadata.mutation_str,
        ] + self.result_energy_rotamer

class EnergyCsvStreams:

    class Stream(NamedTuple):
        stream : TextIO
        writer : '_writer'

    K_POTENTIAL = "potential"
    K_ELECTROSTATICS = "electrostatics"
    K_ELECTROSTATICS_SR = "electrostatics.sr"
    K_LJ = "vanderwaals"
    K_LJ_SR = "vanderwaals.sr"

    STREAMS = [K_POTENTIAL, K_ELECTROSTATICS, K_ELECTROSTATICS_SR, K_LJ, K_LJ_SR]

    def __init__(self, directory: str, suffix : str = "summary"):
        self.__directory = directory
        self.__suffix = suffix
        self.__streams : Dict[str, EnergyCsvStreams.Stream] = {}

    def __file(self, component : str) -> str:
        return path.join(self.__directory, "%s.%s.csv" % (component, self.__suffix))

    def __enter__(self, *args, **kwargs):

        try:
            for stream_name in self.STREAMS:
                stream = open(self.__file(stream_name), 'w').__enter__(*args, **kwargs)
                self.__streams[stream_name] = self.Stream(
                    stream = stream,
                    writer = csv.writer(stream)
                )

            return self
        finally:
            self.__exit__()

    def __exit__(self, *args, **kwargs):

        exception = None
        for stream in self.__streams.values():
            try:
                stream.stream.__exit__(*args, **kwargs)
            except Exception as e:
                exception = e
        
        self.__streams = {}

        if exception:
            raise exception

    def write(self, result : 'EnergyMinimizationCombined'):
        metadata = result.result_metadata
        common_columns : List[Any] = [
            metadata.structure_name,
            metadata.mutation_str
        ]

        rotamers_entry : List[Any] = [r for r in result.result_energy_rotamer if r is not None]

        csv_rows = {
            self.K_POTENTIAL: [r.potential for r in rotamers_entry],
            self.K_ELECTROSTATICS: [r.electrostatic for r in rotamers_entry],
            self.K_ELECTROSTATICS_SR: [r.electrostatic_sr for r in rotamers_entry],
            self.K_LJ: [r.lj for r in rotamers_entry],
            self.K_LJ_SR: [r.lj_sr for r in rotamers_entry]
        }

        for k, csv_row in csv_rows.items():
            row = common_columns + csv_row
            self.__streams[k].writer.writerow(row)

def energy_runs_to_csv(
    directory : str,
    result_directory : Optional[str] = None
):

    result_directory = result_directory or directory
    runs = EnergyMinimizationRun.from_folders(directory)

    with EnergyCsvStreams(result_directory) as csv_energy:
        for result in EnergyMinimizationCombined.summarize_results(runs):
            csv_energy.write(result)