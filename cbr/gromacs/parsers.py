from typing import Iterable, List, NamedTuple, Optional, TextIO

def split_ixs(line : str, indexes : Iterable[int]):
    ix_0 = 0
    for index in indexes:
        yield line[ix_0:index + 1].strip()
        ix_0 = index + 1

def parse_energy_to_json(log : TextIO) -> List[dict]:

    result = []
    current_line = log.readline()

    def stream_done():
        return current_line == ''
    while(not stream_done()):
        if "Energies" in current_line:
            current_line = log.readline()
            values = []

            while not "Energies" in current_line and not stream_done():

                try:
                    headers_line = current_line
                    current_line = entries_line = log.readline()
                    entries = entries_line.split()
                    entries_values = map(float, entries)
                    indexes = [entries_line.find(entry) + len(entry) - 1 for entry in entries]
                    headers = split_ixs(headers_line, indexes)

                    values = values + list(zip(headers, map(float, entries_values)))
                except ValueError:
                    pass

            result.append(dict(values))

        else:
            current_line = log.readline()

    return result

class EnergyLogEntry(NamedTuple):
    step : Optional[int]
    potential : float

    @staticmethod
    def from_dict(entry_dict : dict):

        K_POTENTIAL = "Potential"
        K_STEP = "Step"
        step = entry_dict.get(K_STEP)

        return EnergyLogEntry(
            step = step and int(step),
            potential = entry_dict[K_POTENTIAL]
        )

    @staticmethod
    def parse_entries(stream : TextIO) -> 'List[EnergyLogEntry]':
        return list(map(EnergyLogEntry.from_dict, parse_energy_to_json(stream)))