from typing import Iterable, List, NamedTuple, Optional

from ...core.streams import get_input_stream, InputStream
from ..raspp.contacts import Contacts, read_contacts_objects

ChimeraId = List[int]

class SchemaEnergy(NamedTuple):
    class Entry(NamedTuple):
        chimera : Optional[ChimeraId]
        energy : float
        mutations : float

    average_energy : float
    average_mutations : float
    entries : Iterable[Entry]

def parse_schema_energy(in_data : InputStream) -> SchemaEnergy:

    with get_input_stream(in_data) as input:

        avg_e = 0
        avg_m = 0
        results = []

        for line in input.stream.readlines():
            line = line.strip().split()

            if line[0:4] == ['#', 'Average', 'disruption', '<E>']:
                avg_e = float(line[5])
            elif line[0:4] == ['#', 'Average', 'mutation', '<m>']:
                avg_m = float(line[5])
            elif len(line) == 3:

                results.append(
                    SchemaEnergy.Entry(
                        chimera=[int(i) for i in line[0].strip()],
                        energy=float(line[1]),
                        mutations=float(line[2])
                    )
                )

        return SchemaEnergy(
            average_energy=avg_e,
            average_mutations=avg_m,
            entries=results
        )

def parse_schema_contacts(in_data : InputStream) -> Contacts:

    with get_input_stream(in_data) as input:
        return read_contacts_objects(input.stream)
