import re
from typing import Dict, Iterable, List, NamedTuple, Optional

from ...core.streams import get_input_stream, InputStream
from ..raspp.contacts import Contacts, read_contacts_objects

ChimeraId = List[int]

class SchemaEnergy(NamedTuple):
    class Entry(NamedTuple):
        chimera : Optional[ChimeraId]
        energy : float
        mutations : float
        extra : Dict[str, float]

    average_energy : float
    average_mutations : float
    entries : Iterable[Entry]

AVG_PARSER = re.compile(r'#\s+Average\s+(?P<entry>(disruption\s+\<E\>|mutation\s+\<m\>))\s*=\s*(?P<value>(\d+(\.\d+)?))', re.IGNORECASE)

def create_entry(headers: List[str], values: List[str]) -> SchemaEnergy.Entry:

    chimera = None
    energy = None
    mutations = None
    extra = {}
    for (header, value) in zip(headers, values):

        if header == "# chimera":
            chimera = [int(i) for i in value]
        elif header == "E":
            energy = float(value)
        elif header == "m":
            mutations = float(value)
        else:
            extra[header] = float(value)

    line = " ".join(values)
    assert chimera is not None, f"Line has no chimera: {line}"
    assert energy is not None, f"Line has no disruption {line}"
    assert mutations is not None, f"Line has no mutations {line}"

    return SchemaEnergy.Entry(
        chimera=chimera,
        energy=energy,
        mutations=mutations,
        extra=extra
    )


def parse_schema_energy(in_data : InputStream) -> SchemaEnergy:

    with get_input_stream(in_data) as input:

        avg_e = 0
        avg_m = 0
        results : List[SchemaEnergy.Entry] = []
        headers = None

        for line in input.stream.readlines():
            line = line.strip()

            matches_avg = AVG_PARSER.search(line)

            if matches_avg is not None:
                entry = matches_avg.group('entry')
                value = float(matches_avg.group('value'))

                if entry.startswith('disruption'):
                    avg_e = value
                elif entry.startswith('mutation'):
                    avg_m = value
                else:
                    raise ValueError(f"Unknown entry '{entry}'")
                continue

            line = line.split('\t')

            if len(line) == 0:
                continue

            if headers is None:
                headers = line
                continue

            results.append(create_entry(headers, line))

        return SchemaEnergy(
            average_energy=avg_e,
            average_mutations=avg_m,
            entries=results
        )

def parse_schema_contacts(in_data : InputStream) -> Contacts:

    with get_input_stream(in_data) as input:
        return read_contacts_objects(input.stream)
